/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;

namespace Proteomics.Utilities
{
    public static class Constants
    {
        public const double PROTON_MASS = 1.00727647;
        public const double C12_C13_MASS_DIFFERENCE = 1.0033548378;
        public const double WATER_MONOISOTOPIC_MASS = 18.0105646942;
        public const double WATER_AVERAGE_MASS = 18.01528;
        public const double PEPTIDE_N_TERMINAL_MONOISOTOPIC_MASS = 1.0078250321;
        public const double PEPTIDE_N_TERMINAL_AVERAGE_MASS = 1.00794;
        public const double PEPTIDE_C_TERMINAL_MONOISOTOPIC_MASS = 17.0027396621;
        public const double PEPTIDE_C_TERMINAL_AVERAGE_MASS = 17.00734;
    }
    /// <summary>
    /// Mass Tolerance Units enum, as expressed in the Morpheus code
    /// </summary>
    public enum MassToleranceUnits
    {
        Da,
        ppm
    }

    /// <summary>
    /// Static Routines used in Proteomic studies. 
    /// TODO In .Net 4.5, you can specify InLining, which should be benefical for these methods
    /// </summary>
    public static class Numerics
    {
        public static double StandardDeviation(IEnumerable<double> values)
        {
            return MathNet.Numerics.Statistics.Statistics.StandardDeviation(values);
        }

        public static double Variance(IEnumerable<double> values)
        {
            return MathNet.Numerics.Statistics.Statistics.Variance(values);
        }

        public static double MassFromMZ(double mz, int charge)
        {
            return mz * Math.Abs(charge) - charge * Constants.PROTON_MASS;
        }

        public static double MZFromMass(double mass, int charge)
        {
            return (mass + charge * Constants.PROTON_MASS) / Math.Abs(charge);
        }

        public static double MZFromMzSingleCharge(double mz, int charge)
        {
            return (mz + (charge - 1) * Constants.PROTON_MASS) / Math.Abs(charge);
        }

        public static double IsotopicMassShift(int nbIsotope, int charge)
        {
            return (Constants.C12_C13_MASS_DIFFERENCE * nbIsotope) / (double)charge;
        }

        public static double CalculateMassError(double experimental, double theoretical, MassToleranceUnits massErrorUnits)
        {
            if (massErrorUnits == MassToleranceUnits.Da)
            {
                return experimental - theoretical;
            }
            else if (massErrorUnits == MassToleranceUnits.ppm)
            {
                return (experimental - theoretical) / theoretical * 1e6;
            }
            else
            {
                return double.NaN;
            }
        }
        
        public static double MzDifference(double mz1, double mz2, MassToleranceUnits tol)
        {
            return Math.Abs(CalculateMassError(mz1, mz2, tol));
        }//*/
        /*
        public static double MzTop(double mz1, MassTolerance tol)
        {
            return MassTolerance.MzTop(mz1, tol);
        }//*/
    }
}