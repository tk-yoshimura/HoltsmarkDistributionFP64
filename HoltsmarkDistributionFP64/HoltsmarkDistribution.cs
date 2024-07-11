// Author and Approximation Formula Coefficient Generator: T.Yoshimura
// Github: https://github.com/tk-yoshimura
// Original Code: https://github.com/tk-yoshimura/HoltsmarkDistributionFP64
// Yoshimura, T., "Numerical Evaluation and High Precision Approximation Formula for Holtsmark Distribution", June 2024.
// https://www.researchgate.net/publication/382127111_Numerical_Evaluation_and_High_Precision_Approximation_Formula_for_Holtsmark_Distribution

using HoltsmarkDistributionFP64.InternalUtils;
using HoltsmarkDistributionFP64.RandomGeneration;
using System.Collections.ObjectModel;
using System.Diagnostics;
using static System.Double;

namespace HoltsmarkDistributionFP64 {
    [DebuggerDisplay("{ToString(),nq}")]
    public class HoltsmarkDistribution {

        public double Mu { get; }

        public double C { get; }

        private readonly double c_inv;

        private static readonly double entropy_base = 2.06944850513462440032;

        public HoltsmarkDistribution() : this(mu: 0d, c: 1d) { }

        public HoltsmarkDistribution(double c) : this(mu: 0d, c: c) { }

        public HoltsmarkDistribution(double mu, double c) {
            if (!IsFinite(mu)) {
                throw new ArgumentOutOfRangeException(nameof(mu), "Invalid location parameter.");
            }
            if (!(c > 0 && IsFinite(c))) {
                throw new ArgumentOutOfRangeException(nameof(c), "Invalid scale parameter.");
            }

            Mu = mu;
            C = c;

            c_inv = 1d / c;
        }

        public double PDF(double x) {
            double u = (x - Mu) * c_inv;

            if (IsNaN(u)) {
                return NaN;
            }
            if (IsInfinity(u)) {
                return 0d;
            }

            double pdf = PDFPade.Value(u) * c_inv;

            return pdf;
        }

        public double CDF(double x, Interval interval = Interval.Lower) {
            double u = (x - Mu) * c_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            double cdf = (interval == Interval.Lower) ? CDFPade.Value(-u) : CDFPade.Value(u);

            return cdf;
        }

        public double Quantile(double p, Interval interval = Interval.Lower) {
            if (!(p >= 0d && p <= 1d)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                double x = Mu - C * QuantilePade.Value(p);

                return x;
            }
            else {
                double x = Mu + C * QuantilePade.Value(p);

                return x;
            }
        }

        public double Sample(Random random) {
            double u = random.NextUniformOpenInterval01() - 0.5d;
            double w = random.NextUniformOpenInterval0();

            double cu = CosPi(u);

            double r = SinPi(u * 1.5d) * Cbrt(Log(w) / (CosPi(u * 0.5) * cu * cu));
            double v = r * C + Mu;

            return v;
        }

        public bool Symmetric => true;

        public double Median => Mu;

        public double Mode => Mu;

        public double Mean => Mu;

        public double Variance => PositiveInfinity;

        public double Skewness => NaN;

        public double Kurtosis => NaN;

        public double Entropy => entropy_base + Log(C);

        public double Alpha => 1.5d;

        public double Beta => 0d;

        public static HoltsmarkDistribution operator +(HoltsmarkDistribution dist1, HoltsmarkDistribution dist2) {
            return new(dist1.Mu + dist2.Mu, ExMath.Pow2d3(ExMath.Pow3d2(dist1.C) + ExMath.Pow3d2(dist2.C)));
        }

        public static HoltsmarkDistribution operator -(HoltsmarkDistribution dist1, HoltsmarkDistribution dist2) {
            return new(dist1.Mu - dist2.Mu, ExMath.Pow2d3(ExMath.Pow3d2(dist1.C) + ExMath.Pow3d2(dist2.C)));
        }

        public static HoltsmarkDistribution operator +(HoltsmarkDistribution dist, double s) {
            return new(dist.Mu + s, dist.C);
        }

        public static HoltsmarkDistribution operator -(HoltsmarkDistribution dist, double s) {
            return new(dist.Mu - s, dist.C);
        }

        public static HoltsmarkDistribution operator *(HoltsmarkDistribution dist, double k) {
            return new(dist.Mu * k, dist.C * k);
        }

        public static HoltsmarkDistribution operator /(HoltsmarkDistribution dist, double k) {
            return new(dist.Mu / k, dist.C / k);
        }

        public override string ToString() {
            return $"{typeof(HoltsmarkDistribution).Name}[mu={Mu},c={C}]";
        }

        public string Formula => "p(x; mu, c) := stable_distribution(x; alpha = 3/2, beta = 0, mu, c)";

        private static class PDFPade {
            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_0_1 = new(
                new ReadOnlyCollection<double>([
                    2.87352751452164445024e-1,
                    1.18577398160636011811e-3,
                    -2.16526599226820153260e-2,
                    2.06462093371223113592e-3,
                    2.43382128013710116747e-3,
                    -2.15930711444603559520e-4,
                    -1.04197836740809694657e-4,
                    1.74679078247026597959e-5,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    4.12654472808214997252e-3,
                    2.93891863033354755743e-1,
                    8.70867222155141724171e-3,
                    3.15027515421842640745e-2,
                    2.11141832312672190669e-3,
                    1.23545521355569424975e-3,
                    1.58181113865348637475e-4,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_1_2 = new(
                new ReadOnlyCollection<double>([
                    2.02038159607840130389e-1,
                    -1.20368541260123112191e-2,
                    -3.19235497414059987151e-3,
                    8.88546222140257289852e-3,
                    -5.37287599824602316660e-4,
                    -2.39059149972922243276e-4,
                    9.19551014849109417931e-5,
                    -8.45210544648986348854e-6,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    6.11634701234079515138e-1,
                    4.39922162828115412952e-1,
                    1.73609068791154078128e-1,
                    6.15831808473403962054e-2,
                    1.64364949550314788638e-2,
                    2.94399615562137394932e-3,
                    4.99662797033514776061e-4,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_2_4 = new(
                new ReadOnlyCollection<double>([
                    8.45396231261375200568e-2,
                    -9.15509628797205847643e-3,
                    1.82052933284907579374e-2,
                    -2.44157914076021125182e-4,
                    8.40871885414177705035e-4,
                    7.26592615882060553326e-5,
                    -1.87768359214600016641e-6,
                    1.65716961206268668529e-6,
                    -1.73979640146948858436e-7,
                    7.24351142163396584236e-9,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    8.88099527896838765666e-1,
                    6.53896948546877341992e-1,
                    2.96296982585381844864e-1,
                    1.14107585229341489833e-1,
                    3.08914671331207488189e-2,
                    7.03139384769200902107e-3,
                    1.01201814277918577790e-3,
                    1.12200113270398674535e-4,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_4_8 = new(
                new ReadOnlyCollection<double>([
                    1.36729417918039395222e-2,
                    1.19749117683408419115e-2,
                    6.26780921592414207398e-3,
                    1.84846137440857608948e-3,
                    3.39307829797262466829e-4,
                    2.73606960463362090866e-5,
                    -1.14419838471713498717e-7,
                    1.64552336875610576993e-8,
                    -7.95501797873739398143e-10,
                    2.55422885338760255125e-11,
                    -4.12196487201928768038e-13,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.61334003864149486454e0,
                    1.28348868912975898501e0,
                    6.36594545291321210154e-1,
                    2.11478937436277242988e-1,
                    4.71550897200311391579e-2,
                    6.64679677197059316835e-3,
                    4.93706832858615742810e-4,
                    9.26919465059204396228e-6,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_8_16 = new(
                new ReadOnlyCollection<double>([
                    1.90649774685568282390e-3,
                    7.43708409389806210196e-4,
                    9.53777347766128955847e-5,
                    3.79800193823252979170e-6,
                    2.84836656088572745575e-8,
                    -1.22715411241721187620e-10,
                    8.56789906419220801109e-13,
                    -4.17784858891714869163e-15,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    7.29383849235788831455e-1,
                    2.16287201867831015266e-1,
                    3.28789040872705709070e-2,
                    2.64660789801664804789e-3,
                    1.03662724048874906931e-4,
                    1.47658125632566407978e-6,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_16_32 = new(
                new ReadOnlyCollection<double>([
                    3.07231582988207590928e-4,
                    5.16108848485823513911e-5,
                    3.05776014220862257678e-6,
                    7.64787444325088143218e-8,
                    7.40426355029090813961e-10,
                    1.57451122102115077046e-12,
                    -2.14505675750572782093e-15,
                    5.11204601013038698192e-18,
                    -9.00826023095223871551e-21,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    3.28966789835486457746e-1,
                    4.46981634258601621625e-2,
                    3.22521297380474263906e-3,
                    1.31985203433890010111e-4,
                    3.01507121087942156530e-6,
                    3.47777238523841835495e-8,
                    1.50780503777979189972e-10,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_32_64 = new(
                new ReadOnlyCollection<double>([
                    5.25741312407933720817e-5,
                    2.34425802342454046697e-6,
                    3.30042747965497652847e-8,
                    1.58564820095683252738e-10,
                    1.54070758384735212486e-13,
                    -8.89232435250437247197e-17,
                    8.14099948000080417199e-20,
                    -4.61828164399178360925e-23,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.23544974283127158019e-1,
                    6.01210465184576626802e-3,
                    1.45390926665383063500e-4,
                    1.80594709695117864840e-6,
                    1.06088985542982155880e-8,
                    2.20287881724613104903e-11,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_limit = new(
                new ReadOnlyCollection<double>([
                    2.99206710301074508455e-1,
                    -8.62469397757826072306e-1,
                    1.74661995423629075890e-1,
                    8.75909164947413479137e-1,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    -6.07405848111002255020e0,
                    1.34068401972703571636e1,
                ])
            );

            public static double Value(double x) {
                x = Abs(x);

                double y;
                if (x <= 1d) {
                    Debug.WriteLine("pade minimum segment passed");

                    y = ApproxUtil.Pade(x, pade_plus_0_1);
                }
                else if (x <= 2d) {
                    y = ApproxUtil.Pade(x - 1d, pade_plus_1_2);
                }
                else if (x <= 4d) {
                    y = ApproxUtil.Pade(x - 2d, pade_plus_2_4);
                }
                else if (x <= 8d) {
                    y = ApproxUtil.Pade(x - 4d, pade_plus_4_8);
                }
                else if (x <= 16d) {
                    y = ApproxUtil.Pade(x - 8d, pade_plus_8_16);
                }
                else if (x <= 32d) {
                    y = ApproxUtil.Pade(x - 16d, pade_plus_16_32);
                }
                else if (x <= 64d) {
                    y = ApproxUtil.Pade(x - 32d, pade_plus_32_64);
                }
                else {
                    double u = 1d / ExMath.Pow3d2(x);

                    y = ApproxUtil.Pade(u, pade_plus_limit) * u / x;
                }

                return y;
            }
        }

        private static class CDFPade {
            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_0_0p5 = new(
                new ReadOnlyCollection<double>([
                    5.00000000000000000000e-1,
                    -1.34752580674786639030e-1,
                    1.86318418252163378528e-2,
                    1.04499798132512381447e-2,
                    -1.60831910014592923855e-3,
                    1.38823662364438342844e-4,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    3.05200341554753776087e-1,
                    2.12663999430421346175e-1,
                    7.23836000984872591553e-2,
                    1.67941072412796299986e-2,
                    4.71213644318790580839e-3,
                    5.86825130959777535991e-4,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_0p5_1 = new(
                new ReadOnlyCollection<double>([
                    3.60595773518728397351e-1,
                    5.75238626843218819756e-1,
                    -3.31245319943021227117e-1,
                    1.48132966310216368831e-1,
                    -2.32875122617713403365e-2,
                    2.08038303148835575624e-3,
                    6.01511310581302829460e-6,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    2.32264360456739861886e0,
                    6.39715443864749851087e-1,
                    5.03940458163958921325e-1,
                    8.84780893031413729292e-2,
                    3.01497774031208621961e-2,
                    3.45886005612108195390e-3,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_1_2 = new(
                new ReadOnlyCollection<double>([
                    2.43657975600729535515e-1,
                    -6.02286263626532324632e-2,
                    4.68361231392743283350e-2,
                    -1.13497179885838883972e-3,
                    1.20141595689136205012e-3,
                    3.02402304689333413256e-4,
                    -1.22652173865646814676e-6,
                    2.29521832683440044997e-6,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    5.82002427359748247121e-1,
                    3.96529686558825119743e-1,
                    1.49690294526117385174e-1,
                    5.15049953937764895435e-2,
                    1.30218216530450637564e-2,
                    2.53640337919037463659e-3,
                    3.79575042317720710311e-4,
                    2.94034997185982139717e-5,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_2_4 = new(
                new ReadOnlyCollection<double>([
                    1.05039829654829164883e-1,
                    1.66621813028423002562e-2,
                    2.93820049104275137099e-2,
                    3.36850260303189378587e-3,
                    2.27925819398326978014e-3,
                    1.66394162680543987783e-4,
                    4.51400415642703075050e-5,
                    2.12164734714059446913e-7,
                    1.69306881760242775488e-8,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    9.63461239051296108254e-1,
                    6.54183344973801096611e-1,
                    2.92007762594247903696e-1,
                    1.00918751132022401499e-1,
                    2.55899135910670703945e-2,
                    4.85740416919283630358e-3,
                    6.11435190489589619906e-4,
                    4.10953248859973756440e-5,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_4_8 = new(
                new ReadOnlyCollection<double>([
                    3.05754562114095142887e-2,
                    3.25462617990002726083e-2,
                    1.78205524297204753048e-2,
                    5.61565369088816402420e-3,
                    1.05695297340067353106e-3,
                    9.93588579804511250576e-5,
                    2.94302107205379334662e-6,
                    1.09016076876928010898e-8,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.51164395622515150122e0,
                    1.09391911233213526071e0,
                    4.77950346062744800732e-1,
                    1.34082684956852773925e-1,
                    2.37572579895639589816e-2,
                    2.41806218388337284640e-3,
                    1.10378140456646280084e-4,
                    1.31559373832822136249e-6,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_8_16 = new(
                new ReadOnlyCollection<double>([
                    9.47408470248235718880e-3,
                    4.70888722333356024081e-3,
                    8.66397831692913140221e-4,
                    7.11721056656424862090e-5,
                    2.56320582355149253994e-6,
                    3.37749186035552101702e-8,
                    8.32182844837952178153e-11,
                    -8.80541360484428526226e-14,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    6.98261117346347123707e-1,
                    1.97823959738695249267e-1,
                    2.89311735096848395080e-2,
                    2.30087055379997473849e-3,
                    9.60592522700377510007e-5,
                    1.84474415187428058231e-6,
                    1.14339998084523151203e-8,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_16_32 = new(
                new ReadOnlyCollection<double>([
                    3.19610991747326729867e-3,
                    5.11880074251341162590e-4,
                    2.80704092977662888563e-5,
                    6.31310155466346114729e-7,
                    5.29618446795457166842e-9,
                    9.20292337847562746519e-12,
                    -9.16761719448360345363e-15,
                    1.20433396121606479712e-17,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    2.56283944667056551858e-1,
                    2.56811818304462676948e-2,
                    1.26678062261253559927e-3,
                    3.17001344827541091252e-5,
                    3.68737201224811007437e-7,
                    1.47625352605312785910e-9,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_32_64 = new(
                new ReadOnlyCollection<double>([
                    1.11172037056341397612e-3,
                    7.84545643188695076893e-5,
                    1.94862940242223222641e-6,
                    2.02704958737259525509e-8,
                    7.99772378955335076832e-11,
                    6.62544230949971310060e-14,
                    -3.18234118727325492149e-17,
                    2.03424457039308806437e-20,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.17861198759233241198e-1,
                    5.45962263583663240699e-3,
                    1.25274651876378267111e-4,
                    1.46857544539612002745e-6,
                    8.06441204620771968579e-9,
                    1.53682779460286464073e-11,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_limit = new(
                new ReadOnlyCollection<double>([
                    1.99471140200716338970e-1,
                    -6.90933799347184400422e-1,
                    4.30385245884336871950e-1,
                    3.52790131116013716885e-1,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    -5.05959751628952574534e0,
                    8.04408113719341786819e0,
                ])
            );

            public static double Value(double x) {
                if (IsNegative(x)) {
                    return 1d - Value(-x);
                }

                double y;
                if (x <= 0.5d) {
                    Debug.WriteLine("pade minimum segment passed");

                    y = ApproxUtil.Pade(x, pade_plus_0_0p5);
                }
                else if (x <= 1d) {
                    y = ApproxUtil.Pade(x - 0.5d, pade_plus_0p5_1);
                }
                else if (x <= 2d) {
                    y = ApproxUtil.Pade(x - 1d, pade_plus_1_2);
                }
                else if (x <= 4d) {
                    y = ApproxUtil.Pade(x - 2d, pade_plus_2_4);
                }
                else if (x <= 8d) {
                    y = ApproxUtil.Pade(x - 4d, pade_plus_4_8);
                }
                else if (x <= 16d) {
                    y = ApproxUtil.Pade(x - 8d, pade_plus_8_16);
                }
                else if (x <= 32d) {
                    y = ApproxUtil.Pade(x - 16d, pade_plus_16_32);
                }
                else if (x <= 64d) {
                    y = ApproxUtil.Pade(x - 32d, pade_plus_32_64);
                }
                else {
                    double u = 1d / ExMath.Pow3d2(x);

                    y = ApproxUtil.Pade(u, pade_plus_limit) * u;
                }

                return y;
            }
        }

        private static class QuantilePade {
            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_expm1_2 = new(
                new ReadOnlyCollection<double>([
                    0.00000000000000000000e0,
                    7.59789769759814986929e-1,
                    1.27515008642985381862e0,
                    4.38619247097275579086e-1,
                    -1.25521537863031799276e-1,
                    -2.58555599127223857177e-2,
                    1.20249932437303932411e-2,
                    -1.36753104188136881229e-3,
                    6.57491277860092595148e-5,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    2.48696501912062288766e0,
                    2.06239370128871696850e0,
                    5.67577904795053902651e-1,
                    -2.89022828087034733385e-2,
                    -2.17207943286085236479e-2,
                    3.14098307020814954876e-4,
                    3.51448381406676891012e-4,
                    5.71995514606568751522e-5,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_expm2_3 = new(
                new ReadOnlyCollection<double>([
                    3.84521387984759064238e-1,
                    4.15763727809667641126e-1,
                    -1.73610240124046440578e-2,
                    -3.89915764128788049837e-2,
                    1.07252911248451890192e-2,
                    7.62613727089795367882e-4,
                    -3.11382403581073580481e-4,
                    3.93093062843177374871e-5,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    6.76193897442484823754e-1,
                    3.70953499602257825764e-2,
                    -2.84211795745477605398e-2,
                    2.66146101014551209760e-3,
                    1.85436727973937413751e-3,
                    2.00318687649825430725e-4,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_expm3_4 = new(
                new ReadOnlyCollection<double>([
                    4.46943301497773314460e-1,
                    -1.07267614417424412546e-2,
                    -7.21097021064631831756e-2,
                    2.93948745441334193469e-2,
                    -7.33259305010485915480e-4,
                    -1.38660725579083612045e-3,
                    2.95410432808739478857e-4,
                    -2.88688017391292485867e-5,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    -2.72809429017073648893e-2,
                    -7.85526213469762960803e-2,
                    2.41360900478283465241e-2,
                    3.44597797125179611095e-3,
                    -8.65046428689780375806e-4,
                    -1.04147382037315517658e-4,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_expm4_6 = new(
                new ReadOnlyCollection<double>([
                    4.25344469980677332786e-1,
                    3.42055470008289997369e-2,
                    9.33607217644370441642e-2,
                    4.57057092587794346086e-2,
                    1.16149976708336017542e-2,
                    6.40479797962035786337e-3,
                    1.58526153828271386329e-3,
                    3.84032908993313260466e-4,
                    6.98960839033991110525e-5,
                    9.66690587477825432174e-6,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.60044610004497775009e-1,
                    2.41675490962065446592e-1,
                    1.13752642382290596388e-1,
                    4.05058759031434785584e-2,
                    1.59432816225295660111e-2,
                    4.79286678946992027479e-3,
                    1.16048151070154814260e-3,
                    2.01755520912887201472e-4,
                    2.82884561026909054732e-5,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_expm6_8 = new(
                new ReadOnlyCollection<double>([
                    3.68520435599726877886e-1,
                    8.26682725061327242371e-1,
                    6.85235826889543887309e-1,
                    3.28640408399661746210e-1,
                    9.04801242897407528807e-2,
                    1.57470088502958130451e-2,
                    1.61541023176880542598e-3,
                    9.78919203915954346945e-5,
                    9.71371309261213597491e-8,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    2.29132755303753682133e0,
                    1.95530118226232968288e0,
                    9.55029685883545321419e-1,
                    2.68254036588585643328e-1,
                    4.61398419640231283164e-2,
                    4.66131710581568432246e-3,
                    2.94491397241310968725e-4,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_expm8_16 = new(
                new ReadOnlyCollection<double>([
                    3.48432718168951419458e-1,
                    2.99680703419193973028e-1,
                    1.09531896991852433149e-1,
                    2.28766133215975559897e-2,
                    3.09836969941710802698e-3,
                    2.89346186674853481383e-4,
                    1.96344583080243707169e-5,
                    9.48415601271652569275e-7,
                    3.08821091232356755783e-8,
                    5.58003465656339818416e-10,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    8.73938978582311007855e-1,
                    3.21771888210250878162e-1,
                    6.70432401844821772827e-2,
                    9.05369648218831664411e-3,
                    8.50098390828726795296e-4,
                    5.73568804840571459050e-5,
                    2.78374120155590875053e-6,
                    9.03427646135263412003e-8,
                    1.63556457120944847882e-9,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_expm16_32 = new(
                new ReadOnlyCollection<double>([
                    3.41419813138786920868e-1,
                    1.30219412019722274099e-1,
                    2.36047671342109636195e-2,
                    2.67913051721210953893e-3,
                    2.10896260337301129968e-4,
                    1.19804595761611765179e-5,
                    4.91470756460287578143e-7,
                    1.38299844947707591018e-8,
                    2.25766283556816829070e-10,
                    -8.46510608386806647654e-18,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    3.81461950831351846380e-1,
                    6.91390438866520696447e-2,
                    7.84798596829449138229e-3,
                    6.17735117400536913546e-4,
                    3.50937328177439258136e-5,
                    1.43958654321452532854e-6,
                    4.05109749922716264456e-8,
                    6.61306247924109415113e-10,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_expm32_64 = new(
                new ReadOnlyCollection<double>([
                    3.41392032051575965049e-1,
                    1.53372256183388434238e-1,
                    3.33822240038718319714e-2,
                    4.66328786929735228532e-3,
                    4.67981207864367711082e-4,
                    3.48119463063280710691e-5,
                    2.17755850282052679342e-6,
                    7.40424342670289242177e-8,
                    4.61294046336533026640e-9,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    4.49255524669251621744e-1,
                    9.77826688966262423974e-2,
                    1.36596271675764346980e-2,
                    1.37080296105355418281e-3,
                    1.01970588303201339768e-4,
                    6.37846903580539445994e-6,
                    2.16883897125962281968e-7,
                    1.35121503608967367232e-8,
                ])
            );

            public static double Value(double x) {
                if (x > 0.5) {
                    return -Value(1d - x);
                }

                double v;
                int exponent = ILogB(x);

                if (exponent >= -2) {
                    v = ApproxUtil.Pade(-Log2(ScaleB(x, 1)), pade_plus_expm1_2);
                }
                else if (exponent >= -3) {
                    v = ApproxUtil.Pade(-Log2(ScaleB(x, 2)), pade_plus_expm2_3);
                }
                else if (exponent >= -4) {
                    v = ApproxUtil.Pade(-Log2(ScaleB(x, 3)), pade_plus_expm3_4);
                }
                else if (exponent >= -6) {
                    v = ApproxUtil.Pade(-Log2(ScaleB(x, 4)), pade_plus_expm4_6);
                }
                else if (exponent >= -8) {
                    v = ApproxUtil.Pade(-Log2(ScaleB(x, 6)), pade_plus_expm6_8);
                }
                else if (exponent >= -16) {
                    v = ApproxUtil.Pade(-Log2(ScaleB(x, 8)), pade_plus_expm8_16);
                }
                else if (exponent >= -32) {
                    v = ApproxUtil.Pade(-Log2(ScaleB(x, 16)), pade_plus_expm16_32);
                }
                else if (exponent >= -64) {
                    v = ApproxUtil.Pade(-Log2(ScaleB(x, 32)), pade_plus_expm32_64);
                }
                else {
                    v = 1d / ScaleB(Cbrt(Pi), 1);
                }

                double y = v / ExMath.Pow2d3(x);

                return y;
            }
        }
    }
}