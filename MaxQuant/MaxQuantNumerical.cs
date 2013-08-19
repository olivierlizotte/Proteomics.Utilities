/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;

namespace Proteomics.Utilities.MaxQuant
{
    public delegate void LfitFunc(double x, double[] a);

    public delegate void MrqminFunc(double x, double[] a, out double y, double[] dyda, int na);

    /// <summary>
    /// This is a collection of static routines which are adaptations of algorithms from the book 
    /// 'Numerical Recipes in C'. They are almost literal, with a few exceptions: 1) Indices are 
    /// zero-based, not one-based as in the book. 2) Some integer arguments that indicate sizes 
    /// of arrays were omitted. Their values are infered instead directly from the corresponding arrays.
    /// 3) Names of methods are capitalized. 4) Everything is supposed to be thread-safe. Therefore static
    /// fields for temporary storage of values are avoided. 5) Variable declarations may have been 
    /// moved closer to their usage.
    /// </summary>
    public static class NumericalRecipes
    {
        private const double GcfEps = 3.0e-7;
        private const double GcfFpmin = 1.0e-30;
        private const double GcfItmax = 100;
        private const double GserEps = 3.0e-7;
        private const double GserItmax = 100;

        private static readonly double[] GammlnCof =
            new double[]
				{
					76.18009172947146, -86.50532032941677,
					24.01409824083091, -1.231739572450155,
					0.1208650973866179e-2, -0.5395239384953e-5
				};

        private const int betacfMAXIT = 100;
        private const double betacfEPS = 3.0e-7;
        private const double betacfFPMIN = 1.0e-30;

        public static double Betacf(double a, double b, double x)
        {
            int m;
            double qab = a + b;
            double qap = a + 1.0;
            double qam = a - 1.0;
            double c = 1.0;
            double d = 1.0 - qab * x / qap;
            if (Math.Abs(d) < betacfFPMIN)
            {
                d = betacfFPMIN;
            }
            d = 1.0 / d;
            double h = d;
            for (m = 1; m <= betacfMAXIT; m++)
            {
                int m2 = 2 * m;
                double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
                d = 1.0 + aa * d;
                if (Math.Abs(d) < betacfFPMIN)
                {
                    d = betacfFPMIN;
                }
                c = 1.0 + aa / c;
                if (Math.Abs(c) < betacfFPMIN)
                {
                    c = betacfFPMIN;
                }
                d = 1.0 / d;
                h *= d * c;
                aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
                d = 1.0 + aa * d;
                if (Math.Abs(d) < betacfFPMIN)
                {
                    d = betacfFPMIN;
                }
                c = 1.0 + aa / c;
                if (Math.Abs(c) < betacfFPMIN)
                {
                    c = betacfFPMIN;
                }
                d = 1.0 / d;
                double del = d * c;
                h *= del;
                if (Math.Abs(del - 1.0) < betacfEPS)
                {
                    break;
                }
            }
            if (m > betacfMAXIT)
            {
                throw new Exception("a or b too big, or MAXIT too small in betacf");
            }
            return h;
        }

        public static double Betai(double a, double b, double x)
        {
            double bt;
            if (x < 0.0 || x > 1.0)
            {
                throw new Exception("Bad x in routine betai");
            }
            if (x == 0.0 || x == 1.0)
            {
                bt = 0.0;
            }
            else
            {
                bt = Math.Exp(Gammln(a + b) - Gammln(a) - Gammln(b) + a * Math.Log(x) + b * Math.Log(1.0 - x));
            }
            if (x < (a + 1.0) / (a + b + 2.0))
            {
                return bt * Betacf(a, b, x) / a;
            }
            return 1.0 - bt * Betacf(b, a, 1.0 - x) / b;
        }

        public static void Covsrt(double[,] covar)
        {
            int ma = covar.GetLength(0);
            int mfit = ma;
            for (int i = mfit; i < ma; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    covar[i, j] = 0.0;
                    covar[j, i] = 0.0;
                }
            }
            int k = mfit - 1;
            for (int j = ma - 1; j >= 0; j--)
            {
                double swap;
                for (int i = 0; i < ma; i++)
                {
                    swap = covar[i, k];
                    covar[i, k] = covar[i, j];
                    covar[i, j] = swap;
                }
                for (int i = 0; i < ma; i++)
                {
                    swap = covar[k, i];
                    covar[k, i] = covar[j, i];
                    covar[j, i] = swap;
                }
                k--;
            }
        }

        /// <summary>
        /// Returns the error function erf(x)
        /// </summary>
        public static double Erff(double x)
        {
            return x < 0.0 ? -Gammp(0.5, x * x) : Gammp(0.5, x * x);
        }

        /// <summary>
        /// Returns the complementary error function erfc(x) = 1 - erf(x)
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double Erffc(double x)
        {
            return x < 0.0 ? 1.0 + Gammp(0.5, x * x) : Gammq(0.5, x * x);
        }

        private static readonly double[] aaa = new double[501];

        public static double Factln(int n)
        {
            if (n < 0)
            {
                throw new Exception("Negative factorial in routine factln");
            }
            if (n <= 1)
            {
                return 0.0;
            }
            if (n <= 500)
            {
                return aaa[n] != 0 ? aaa[n] : (aaa[n] = Gammln(n + 1.0));
            }
            return Gammln(n + 1.0);
        }


        public static double Gammln(double xx)
        {
            double x;
            int j;
            double y = x = xx;
            double tmp = x + 5.5;
            tmp -= (x + 0.5) * Math.Log(tmp);
            double ser = 1.000000000190015;
            for (j = 0; j <= 5; j++)
            {
                ser += GammlnCof[j] / ++y;
            }
            return -tmp + Math.Log(2.5066282746310005 * ser / x);
        }

        /// <summary>
        /// Returns the incomplete gamma function P(a,x)
        /// </summary>
        public static double Gammp(double a, double x)
        {
            double gamser = Double.NaN, gammcf = Double.NaN, gln = Double.NaN;
            if (x < 0.0 || a <= 0.0)
            {
                throw new Exception("Invalid arguments in routine gammq");
            }
            if (x < (a + 1.0))
            {
                Gser(ref gamser, a, x, ref gln);
                return gamser;
            }
            Gcf(ref gammcf, a, x, ref gln);
            return 1.0 - gammcf;
        }

        /// <summary>
        /// Returns the incomplete gamma function Q(a,x) = 1 - P(a,x)
        /// </summary>
        public static double Gammq(double a, double x)
        {
            if (x < 0.0 || a <= 0.0)
            {
                throw new Exception("Invalid arguments in routine gammq");
            }
            if (x < (a + 1.0))
            {
                double gamser = Double.NaN, gln = Double.NaN;
                Gser(ref gamser, a, x, ref gln);
                return 1.0 - gamser;
            }
            else
            {
                double gammcf = Double.NaN, gln = Double.NaN;
                Gcf(ref gammcf, a, x, ref gln);
                return gammcf;
            }
        }

        public static void Gaussj(double[,] a, int n, double[,] b, int m)
        {
            int i, icol = -1, irow = -1, j, k, l;
            double temp;
            int[] indxc = new int[n];
            int[] indxr = new int[n];
            int[] ipiv = new int[n];
            for (j = 0; j < n; j++)
            {
                ipiv[j] = 0;
            }
            for (i = 0; i < n; i++)
            {
                double big = 0.0;
                for (j = 0; j < n; j++)
                {
                    if (ipiv[j] != 1)
                    {
                        for (k = 0; k < n; k++)
                        {
                            if (ipiv[k] == 0)
                            {
                                if (Math.Abs(a[j, k]) >= big)
                                {
                                    big = Math.Abs(a[j, k]);
                                    irow = j;
                                    icol = k;
                                }
                            }
                            else if (ipiv[k] > 1)
                            {
                                throw new Exception("gaussj: Singular Matrix-1");
                            }
                        }
                    }
                }
                ++(ipiv[icol]);
                if (irow != icol)
                {
                    for (l = 0; l < n; l++)
                    {
                        temp = a[irow, l];
                        a[irow, l] = a[icol, l];
                        a[icol, l] = temp;
                    }
                    for (l = 0; l < m; l++)
                    {
                        temp = b[irow, l];
                        b[irow, l] = b[icol, l];
                        b[icol, l] = temp;
                    }
                }
                indxr[i] = irow;
                indxc[i] = icol;
                if (a[icol, icol] == 0.0)
                {
                    throw new Exception("gaussj: Singular Matrix-2");
                }
                double pivinv = 1.0 / a[icol, icol];
                a[icol, icol] = 1.0;
                for (l = 0; l < n; l++)
                {
                    a[icol, l] *= pivinv;
                }
                for (l = 0; l < m; l++)
                {
                    b[icol, l] *= pivinv;
                }
                int ll;
                for (ll = 0; ll < n; ll++)
                {
                    if (ll != icol)
                    {
                        double dum = a[ll, icol];
                        a[ll, icol] = 0.0;
                        for (l = 0; l < n; l++)
                        {
                            a[ll, l] -= a[icol, l] * dum;
                        }
                        for (l = 0; l < m; l++)
                        {
                            b[ll, l] -= b[icol, l] * dum;
                        }
                    }
                }
            }
            for (l = n - 1; l >= 0; l--)
            {
                if (indxr[l] != indxc[l])
                {
                    for (k = 0; k < n; k++)
                    {
                        temp = a[k, indxr[l]];
                        a[k, indxr[l]] = a[k, indxc[l]];
                        a[k, indxc[l]] = temp;
                    }
                }
            }
        }

        public static void Gcf(ref double gammcf, double a, double x, ref double gln)
        {
            int i;
            gln = Gammln(a);
            double b = x + 1.0 - a;
            double c = 1.0 / GcfFpmin;
            double d = 1.0 / b;
            double h = d;
            for (i = 1; i <= GcfItmax; i++)
            {
                double an = -i * (i - a);
                b += 2.0;
                d = an * d + b;
                if (Math.Abs(d) < GcfFpmin)
                {
                    d = GcfFpmin;
                }
                c = b + an / c;
                if (Math.Abs(c) < GcfFpmin)
                {
                    c = GcfFpmin;
                }
                d = 1.0 / d;
                double del = d * c;
                h *= del;
                if (Math.Abs(del - 1.0) < GcfEps)
                {
                    break;
                }
            }
            if (i > GcfItmax)
            {
                throw new Exception("a too large, ITMAX too small in gcf");
            }
            gammcf = Math.Exp(-x + a * Math.Log(x) - gln) * h;
        }

        public static void Gser(ref double gamser, double a, double x, ref double gln)
        {
            gln = Gammln(a);
            if (x <= 0.0)
            {
                if (x < 0.0)
                {
                    throw new Exception("x less than 0 in routine gser");
                }
                gamser = 0.0;
                return;
            }
            double ap = a;
            double sum = 1.0 / a;
            double del = sum;
            for (int n = 1; n <= GserItmax; n++)
            {
                ++ap;
                del *= x / ap;
                sum += del;
                if (Math.Abs(del) < Math.Abs(sum) * GserEps)
                {
                    gamser = sum * Math.Exp(-x + a * Math.Log(x) - (gln));
                    return;
                }
            }
            throw new Exception("a too large, ITMAX too small in routine gser");
        }

        public static void Lfit(double[] x, double[] y, double[] sig, double[] a, out double[,] covar, out double chisq,
                                LfitFunc funcs)
        {
            int ndat = x.Length;
            int i, j;
            int l;
            int ma = a.Length;
            int mfit = a.Length;
            double[,] beta = new double[ma, 1];
            double[] afunc = new double[ma];
            covar = new double[ma, ma];
            for (i = 0; i < ndat; i++)
            {
                funcs(x[i], afunc);
                double ym = y[i];
                double sig2i = 1.0 / (sig[i] * sig[i]);
                for (l = 0; l < ma; l++)
                {
                    double wt = afunc[l] * sig2i;
                    for (int m = 0; m <= l; m++)
                    {
                        covar[l, m] += wt * afunc[m];
                    }
                    beta[l, 0] += ym * wt;
                }
            }
            for (j = 1; j < mfit; j++)
            {
                int k;
                for (k = 0; k < j; k++)
                {
                    covar[k, j] = covar[j, k];
                }
            }
            Gaussj(covar, mfit, beta, 1);
            for (l = 0; l < ma; l++)
            {
                a[l] = beta[l, 0];
            }
            chisq = 0.0;
            for (i = 0; i < ndat; i++)
            {
                funcs(x[i], afunc);
                double sum = 0.0;
                for (j = 0; j < ma; j++)
                {
                    sum += a[j] * afunc[j];
                }
                chisq += ((y[i] - sum) / sig[i]) * ((y[i] - sum) / sig[i]);
            }
            Covsrt(covar);
        }

        public static void Mrqcof(double[] x, double[] y, double[] sig, int ndata,
                                  double[] a, double[,] alpha, double[] beta,
                                  out double chisq, MrqminFunc func)
        {
            int ma = a.Length;
            double[] dyda = new double[ma];
            for (int j = 0; j < ma; j++)
            {
                for (int k = 0; k <= j; k++)
                {
                    alpha[j, k] = 0.0;
                }
                beta[j] = 0.0;
            }
            chisq = 0.0;
            for (int i = 0; i < ndata; i++)
            {
                double ymod;
                func(x[i], a, out ymod, dyda, ma);
                double sig2i = 1.0 / (sig[i] * sig[i]);
                double dy = y[i] - ymod;
                for (int j = 0, l = 0; l < ma; l++)
                {
                    double wt = dyda[l] * sig2i;
                    for (int k = 0, m = 0; m <= l; m++)
                    {
                        alpha[j, k++] += wt * dyda[m];
                    }
                    beta[j] += dy * wt;
                    j++;
                }
                chisq += dy * dy * sig2i;
            }
            for (int j = 1; j < ma; j++)
            {
                for (int k = 0; k < j; k++)
                {
                    alpha[k, j] = alpha[j, k];
                }
            }
        }

        public static void Mrqmin(double[] x, double[] y, double[] sig, int ndata, double[] a,
                                  int[] ia, int ma, double[,] covar, double[,] alpha, out double chisq, MrqminFunc func,
                                  ref double alamda, ref double ochisq, ref double[,] oneda, ref int mfit, ref double[] atry,
                                  ref double[] beta, ref double[] da)
        {
            int j, k, l, m;
            if (alamda < 0.0)
            {
                atry = new double[ma];
                beta = new double[ma];
                da = new double[ma];
                mfit = ma;
                oneda = new double[mfit, 1];
                alamda = 0.001;
                Mrqcof(x, y, sig, ndata, a, alpha, beta, out chisq, func);
                ochisq = (chisq);
                for (j = 0; j < ma; j++)
                {
                    atry[j] = a[j];
                }
            }
            for (j = 0, l = 0; l < ma; l++)
            {
                for (k = 0, m = 0; m < ma; m++)
                {
                    covar[j, k] = alpha[j, k];
                    k++;
                }
                covar[j, j] = alpha[j, j] * (1.0 + (alamda));
                oneda[j, 0] = beta[j];
                j++;
            }
            Gaussj(covar, mfit, oneda, 1);
            for (j = 0; j < mfit; j++)
            {
                da[j] = oneda[j, 0];
            }
            if (alamda == 0.0)
            {
                Covsrt(covar);
                chisq = ochisq;
                return;
            }
            for (j = 0, l = 0; l < ma; l++)
            {
                atry[l] = a[l] + da[j++];
            }
            Mrqcof(x, y, sig, ndata, atry, covar, da, out chisq, func);
            if (chisq < ochisq)
            {
                alamda *= 0.1;
                ochisq = (chisq);
                for (j = 0, l = 0; l < ma; l++)
                {
                    for (k = 0, m = 0; m < ma; m++)
                    {
                        alpha[j, k] = covar[j, k];
                        k++;
                    }
                    beta[j] = da[j];
                    a[l] = atry[l];
                    j++;
                }
            }
            else
            {
                alamda *= 10.0;
                chisq = ochisq;
            }
        }

        /// <summary>
        /// Computes (a^2+b^2)^1/2 without destructive underflow or overflow.
        /// </summary>
        public static double Pythag(double a, double b)
        {
            double absa = Math.Abs(a);
            double absb = Math.Abs(b);
            if (absa > absb)
            {
                return absa * Math.Sqrt(1.0 + (absb / absa) * (absb / absa));
            }
            return (absb == 0.0 ? 0.0 : absb * Math.Sqrt(1.0 + (absa / absb) * (absa / absb)));
        }

        public static void Tqli(double[] d, double[] e, double[,] z)
        {
            int n = d.Length;
            int l;
            int i;
            for (i = 1; i < n; i++)
            {
                e[i - 1] = e[i];
            }
            e[n - 1] = 0.0;
            for (l = 0; l < n; l++)
            {
                int iter = 0;
                int m;
                do
                {
                    for (m = l; m < n - 1; m++)
                    {
                        double dd = Math.Abs(d[m]) + Math.Abs(d[m + 1]);
                        if ((Math.Abs(e[m]) + dd) == dd)
                        {
                            break;
                        }
                    }
                    if (m != l)
                    {
                        if (iter++ == 30)
                        {
                            throw new Exception("Too many iterations in tqli");
                        }
                        double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                        double r = Pythag(g, 1.0);
                        g = d[m] - d[l] + e[l] / (g + ((g) >= 0.0 ? Math.Abs(r) : -Math.Abs(r)));
                        double c;
                        double s = c = 1.0;
                        double p = 0.0;
                        for (i = m - 1; i >= l; i--)
                        {
                            double f = s * e[i];
                            double b = c * e[i];
                            e[i + 1] = (r = Pythag(f, g));
                            if (r == 0.0)
                            {
                                d[i + 1] -= p;
                                e[m] = 0.0;
                                break;
                            }
                            s = f / r;
                            c = g / r;
                            g = d[i + 1] - p;
                            r = (d[i] - g) * s + 2.0 * c * b;
                            d[i + 1] = g + (p = s * r);
                            g = c * r - b;
                            int k;
                            for (k = 0; k < n; k++)
                            {
                                f = z[k, i + 1];
                                z[k, i + 1] = s * z[k, i] + c * f;
                                z[k, i] = c * z[k, i] - s * f;
                            }
                        }
                        if (r == 0.0 && i >= l)
                        {
                            continue;
                        }
                        d[l] -= p;
                        e[l] = g;
                        e[m] = 0.0;
                    }
                } while (m != l);
            }
        }

        public static void Tred2(double[,] a, out double[] d, out double[] e)
        {
            int n = a.GetLength(0);
            d = new double[n];
            e = new double[n];
            int l, k, j, i;
            double g;
            for (i = n - 1; i >= 1; i--)
            {
                l = i - 1;
                double scale;
                double h = scale = 0.0;
                if (l > 0)
                {
                    for (k = 0; k <= l; k++)
                    {
                        scale += Math.Abs(a[i, k]);
                    }
                    if (scale == 0.0)
                    {
                        e[i] = a[i, l];
                    }
                    else
                    {
                        for (k = 0; k <= l; k++)
                        {
                            a[i, k] /= scale;
                            h += a[i, k] * a[i, k];
                        }
                        double f = a[i, l];
                        g = (f >= 0.0 ? -Math.Sqrt(h) : Math.Sqrt(h));
                        e[i] = scale * g;
                        h -= f * g;
                        a[i, l] = f - g;
                        f = 0.0;
                        for (j = 0; j <= l; j++)
                        {
                            a[j, i] = a[i, j] / h;
                            g = 0.0;
                            for (k = 0; k <= j; k++)
                            {
                                g += a[j, k] * a[i, k];
                            }
                            for (k = j + 1; k <= l; k++)
                            {
                                g += a[k, j] * a[i, k];
                            }
                            e[j] = g / h;
                            f += e[j] * a[i, j];
                        }
                        double hh = f / (h + h);
                        for (j = 0; j <= l; j++)
                        {
                            f = a[i, j];
                            e[j] = g = e[j] - hh * f;
                            for (k = 0; k <= j; k++)
                            {
                                a[j, k] -= (f * e[k] + g * a[i, k]);
                            }
                        }
                    }
                }
                else
                {
                    e[i] = a[i, l];
                }
                d[i] = h;
            }
            d[0] = 0.0;
            e[0] = 0.0;
            /* Contents of this loop can be omitted if eigenvectors not
            wanted except for statement d[i]=a[i][i]; */
            for (i = 0; i < n; i++)
            {
                l = i - 1;
                if (d[i] != 0)
                {
                    for (j = 0; j <= l; j++)
                    {
                        g = 0.0;
                        for (k = 0; k <= l; k++)
                        {
                            g += a[i, k] * a[k, j];
                        }
                        for (k = 0; k <= l; k++)
                        {
                            a[k, j] -= g * a[k, i];
                        }
                    }
                }
                d[i] = a[i, i];
                a[i, i] = 1.0;
                for (j = 0; j <= l; j++)
                {
                    a[j, i] = a[i, j] = 0.0;
                }
            }
        }
        public static void LinFit2(double[] x, double[] y, double[] a, LfitFunc funcs)
        {
            double chisq;
            LinFit2(x, y, null, a, out chisq, funcs);
        }

        public static void LinFit2(double[] x, double[] y, double[] sig, double[] a,
                                   out double chisq, LfitFunc funcs)
        {
            double[,] covar;
            if (sig == null)
            {
                sig = new double[x.Length];
                for (int i = 0; i < sig.Length; i++)
                {
                    sig[i] = 1E-2;
                }
            }
            Lfit(x, y, sig, a, out covar, out chisq, funcs);
        }

        private const double mass0 = 445.120025;
        public static double RelativeCorrection(double u, double[,] mzCalibrationPar, double[,] intensityCalibrationPar,
                                                double intens, int massRange)
        {
            if (mzCalibrationPar == null)
            {
                return 0;
            }
            u = 100 / Math.Sqrt(u) - 100 / Math.Sqrt(mass0);
            double w = Math.Log(intens) - Math.Log(1e7);
            double result = 0;
            double fact = u;
            for (int i = 0; i < mzCalibrationPar.GetLength(1); i++)
            {
                result += fact * mzCalibrationPar[massRange, i];
                fact *= u;
            }
            fact = w;
            for (int i = 0; i < intensityCalibrationPar.GetLength(1); i++)
            {
                result += fact * intensityCalibrationPar[massRange, i];
                fact *= w;
            }
            return result;
        }
    }
}