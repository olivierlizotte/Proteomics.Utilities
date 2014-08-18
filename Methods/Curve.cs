using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteomics.Utilities.Methods
{
    public class MaxFlowElutionCurve : ElutionCurve
    {
        public int nbProducts;
        public MaxFlowElutionCurve(int nbProductsUsed)
        {
            this.nbProducts = nbProductsUsed;
        }
    }

    public class ElutionCurveMerger
    {
        private List<ElutionCurve> Curves = new List<ElutionCurve>();
        private List<double>       Factor = new List<double>();
        public void AddCurve(ElutionCurve newCurve, double weight)
        {
            Curves.Add(newCurve);
            Factor.Add(weight);
        }

        public ElutionCurve Merge()
        {
            if(Curves.Count > 1)
            {
                ElutionCurve newCurve = new ElutionCurve();
                double sum = 0.0;
                foreach (double val in Factor)
                    sum += val;

                Dictionary<double, int> times = new Dictionary<double, int>();
                foreach (ElutionCurve curve in Curves)
                    foreach (double timePoint in curve.time)
                        if (!times.ContainsKey(timePoint))
                            times.Add(timePoint, 1);
                        else
                            times[timePoint]++;
                List<double> sortedTime = new List<double>(times.Keys);
                sortedTime.Sort();
                foreach (double key in sortedTime)
                    if (times[key] > 1)
                    {
                        double cumulIntensity = 0.0;
                        for (int i = 0; i < Curves.Count; i++)
                            cumulIntensity += Curves[i].InterpolateIntensity(key) * Factor[i] / sum;

                        newCurve.AddPoint(key, cumulIntensity);
                    }
                return newCurve;
            }
            else if (Curves.Count == 1)
                return Curves[0];
            return new ElutionCurve();
        }
    }

    public class ElutionCurve
    {
        public double Area = 0.0;
        public double[] Coefficients = null;        

        public List<double> time = null;
        public List<double> intensityCount = null;

        public double[] GetTimePoints(int nbTimePoints)
        {
            double highestIntensity = 0;
            int indexHighestIntensity = 0;
            for(int i = 0; i < time.Count; i++)
                if(intensityCount[i] > highestIntensity)
                {
                    highestIntensity = intensityCount[i];
                    indexHighestIntensity = i;
                }
                                
            int minIndex = indexHighestIntensity;
            while(minIndex > 0 && !(intensityCount[minIndex] == 0 && intensityCount[minIndex-1] == 0))
                minIndex--;

            int maxIndex = indexHighestIntensity;
            while(maxIndex < time.Count - 1 && !(intensityCount[maxIndex] == 0 && intensityCount[maxIndex + 1] == 0))
                maxIndex++;

            if (maxIndex - minIndex > time.Count)
            {
                double minTime = time[minIndex];
                double maxTime = time[maxIndex];
                double[] points = new double[nbTimePoints];
                for (int i = 0; i < nbTimePoints; i++)
                    points[i] = minTime + (i / (double)nbTimePoints) * (maxTime - minTime);
                return points;
            }
            else
                return time.ToArray();
        }

        public int GetNbPoints()
        {
            return intensityCount.Count;
        }

        private MathNet.Numerics.Interpolation.IInterpolation interpole = null;
        
        public static ElutionCurve Create(Dictionary<double, double> dicOfTimeInMsVsIntensityPerMs)
        {
            ElutionCurve theCurve = new ElutionCurve();
            // -- Test curve fitting function -- //
            theCurve.time = new List<double>();
            theCurve.intensityCount = new List<double>();
            List<double> sortedKeys = new List<double>(dicOfTimeInMsVsIntensityPerMs.Keys);
            sortedKeys.Sort();
            foreach (double time in sortedKeys)
            {
                theCurve.time.Add(time);
                theCurve.intensityCount.Add(dicOfTimeInMsVsIntensityPerMs[time]);
            }
            theCurve.Compute();
            return theCurve;
        }

        public double InterpolateIntensity(double timePoint)
        {
            if (time != null && time.Count > 4)
            {
                if (interpole == null)
                {
                    try
                    {
                        //interpole = new MathNet.Numerics.Interpolation.Algorithms.CubicSplineInterpolation(time, intensityCount);
                        interpole = new MathNet.Numerics.Interpolation.Algorithms.AkimaSplineInterpolation(time, intensityCount);
                        //interpole = MathNet.Numerics.Interpolation.Interpolate.LinearBetweenPoints(time, intensityCount);
                    }catch(Exception ex)
                    {
                        Console.WriteLine(ex.StackTrace);
                    }
                }
                double intensity = interpole.Interpolate(timePoint);
                if (intensity < 0)
                    return 0;
                else
                    return intensity;
            }
            return 0;
        }

        public void Compute()
        {
            if (time != null && time.Count > 4)
            {/*
                List<double> tmpTime = new List<double>(time);
                List<double> tmpIntC = new List<double>(intensityCount);
                
                int nbIter = 2;
                while(nbIter > 0)
                {
                    double worst = 0.0;
                    Area = Proteomics.Utilities.CurveFitter.FitToPolynomial(tmpTime.ToArray(), tmpIntC.ToArray(), out Coefficients);
                    
                }
                //*/
                
                double area1 = Proteomics.Utilities.CurveFitter.FitToPolynomial(time.ToArray(), intensityCount.ToArray(), out Coefficients);
                //double areaInterpol = Proteomics.Utilities.CurveFitter.AreaUnderTheCurve(time, intensityCount);
                double area2 = Proteomics.Utilities.CurveFitter.AreaUnderTheCurve(time, intensityCount);
                //if(area1 / area2 > 2 || area1 / area2 < 0.5)
                //    Console.WriteLine("Too much diff");
                Area = area2;// (area1 + area2) * 0.5;
                //if (Area > areaInterpol * 1.2 || Area < areaInterpol * 0.8)
                //    Console.WriteLine("Too much diff");//*/
            }
            else
            {
                Area = 0;
                Coefficients = new double[0];
            }
        }

        public void AddPoint(double newTimePoint, double newIntensityPerMilliSeconds)
        {
            if (time == null || intensityCount == null)
            {
                time = new List<double>();
                intensityCount = new List<double>();
                time.Add(newTimePoint);
                intensityCount.Add(newIntensityPerMilliSeconds);
            }
            else
                if (time[time.Count - 1] < newTimePoint)
                {
                    time.Add(newTimePoint);
                    intensityCount.Add(newIntensityPerMilliSeconds);
                }
                else
                    Console.WriteLine("UnSorted timepoint inserted");
        }
    }
}
