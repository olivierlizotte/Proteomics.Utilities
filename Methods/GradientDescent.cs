using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteomics.Utilities.Methods
{
    public static class GradientDescent
    {
        public static double ComputeOver(Dictionary<double, double> virtualMixed, Dictionary<double, double> mixed)        
        {
            double cumulOver = 0;
            foreach (double key in virtualMixed.Keys)
                if (virtualMixed[key] > mixed[key])
                    cumulOver += virtualMixed[key] - mixed[key];
            return cumulOver;
        }

        public static double ComputeUnder(Dictionary<double, double> virtualMixed, Dictionary<double, double> mixed)
        {
            double cumulUnder = 0;
            foreach (double key in virtualMixed.Keys)
                if (mixed[key] > virtualMixed[key])
                    cumulUnder += mixed[key] - virtualMixed[key];
            return cumulUnder;
        }
        
        public static void SolveMinimaStyle(List<Dictionary<double, double>> units, Dictionary<double, double> mixed,
                                 out List<double> solution, out double underflow, IConSol ConSole)
        {
            List<double> localFlows = new List<double>();
            foreach (Dictionary<double, double> unit in units)
                localFlows.Add(FindLocalMaxima(unit, mixed));

            Dictionary<double, double> virtualMixed = BuildVirtualDic(localFlows, units, mixed.Count);
            double overError = ComputeOver(virtualMixed, mixed);
            double underError = ComputeUnder(virtualMixed, mixed);

            int bestUnit = 0;
            while (overError >= 1 && bestUnit >= 0)
            {
                //double bestFlow = 0;
                double bestMinima = 0;
                bestUnit = -1;
                for (int i = 0; i < units.Count; i++)
                {
                    if (localFlows[i] > 0)
                    {
                        double minima = FindLocalMinima(units[i], mixed, virtualMixed);
                        //double currentFlow = localFlows[i];
                        //localFlows[i] = minima;
                        //Dictionary<double, double> tmpDic = BuildVirtualDic(localFlows, units);
                        //double tmpUnderError = ComputeUnder(tmpDic, mixed);
                        //double tmpOverError  = ComputeOver(tmpDic, mixed);

                        //double tmpFlowRate = Math.Abs(overError - tmpOverError);
                        //if (tmpUnderError > underError)
                        //    tmpFlowRate /= tmpUnderError - underError;
                        //if (tmpFlowRate > bestFlow)
                        if(minima > bestMinima)
                        {
                            //bestFlow = tmpFlowRate;
                            bestMinima = minima;
                            bestUnit = i;
                        }
                        //localFlows[i] = currentFlow;
                    }
                }
                if (bestUnit >= 0)
                {
                    if (bestMinima > 1)
                        localFlows[bestUnit] -= 1.0;// *0.01
                    else 
                        localFlows[bestUnit] -= bestMinima;// *0.01;
                    if (localFlows[bestUnit] < 0)
                        localFlows[bestUnit] = 0.0;

                    virtualMixed = BuildVirtualDic(localFlows, units, mixed.Count);
                    overError = ComputeOver(virtualMixed, mixed);
                    underError = ComputeUnder(virtualMixed, mixed);
                }
            }//End of while overflow > 1

            solution = new List<double>();
            foreach (double localFlow in localFlows)
                if (overError <= 1.0)
                    solution.Add(localFlow);
                else
                    solution.Add(0);

            underflow = underError;
        }//*/

        public static void SolveMaxFlowStyle(List<Dictionary<double, double>> units, Dictionary<double, double> mixed,
                                 out List<double> solution, out double underflow, IConSol ConSole)
        {
            List<double> localFlows = new List<double>();
            foreach (Dictionary<double, double> unit in units)
                localFlows.Add(FindLocalMaxima(unit, mixed));

            Dictionary<double, double> virtualMixed = BuildVirtualDic(localFlows, units, mixed.Count);
            double overError = ComputeOver(virtualMixed, mixed);
            double underError = ComputeUnder(virtualMixed, mixed);
            double[] bestIndexes = new double[units.Count];

            int iterSize = 1;
            double bestOverallError = double.MaxValue;
            List<double> bestLocalFlows = new List<double>();
            Random rnd = new Random();
            while (overError >= 1 && iterSize < 10000)//anything less than 1 is an acceptable solution
            {
                for (int index = 0; index < bestIndexes.Length; index++)
                    bestIndexes[index] = -1;

                for (int i = 0; i < units.Count; i++)
                {
                    if (localFlows[i] > 0)
                    {
                        localFlows[i] -= iterSize;

                        virtualMixed = BuildVirtualDic(localFlows, units, mixed.Count);
                        double tmpErrorOver = ComputeOver(virtualMixed, mixed);
                        double tmpErrorUnder = ComputeUnder(virtualMixed, mixed);

                        double tmpFlowRate = Math.Abs(overError - tmpErrorOver);
                        double underDiff = 0;
                        if (tmpErrorUnder > underError)
                            underDiff = tmpErrorUnder - underError;
                        if (underDiff >= 1)
                            tmpFlowRate /= underDiff;
                        bestIndexes[i] = tmpFlowRate;

                        localFlows[i] += iterSize;
                    }
                }

                //Pick pseudo randomly best index
                double worstFlowRate = 0.0;
                for (int index = 0; index < bestIndexes.Length; index++)
                    if (bestIndexes[index] > worstFlowRate)
                    {
                        worstFlowRate = bestIndexes[index];
                    }

                if (worstFlowRate > 0)
                {
                    int nbMatching = 0;
                    for (int index = 0; index < bestIndexes.Length; index++)
                        if (bestIndexes[index] >= worstFlowRate)
                            nbMatching++;

                    int iterChoice = rnd.Next(0, nbMatching - 1);
                    int iterNb = 0;
                    for (int index = 0; index < bestIndexes.Length; index++)
                        if (bestIndexes[index] >= worstFlowRate)
                        {
                            if (iterChoice == iterNb)
                            {
                                localFlows[index] -= iterSize;
                                if (localFlows[index] < 0)
                                    localFlows[index] = 0.0;
                            }
                            iterNb++;
                        }
                    iterSize = 1;
                }
                else
                    iterSize++;

                virtualMixed = BuildVirtualDic(localFlows, units, mixed.Count);
                overError = ComputeOver(virtualMixed, mixed);
                underError = ComputeUnder(virtualMixed, mixed);
                if (overError + underError < bestOverallError)
                {
                    bestLocalFlows = new List<double>(localFlows);
                    bestOverallError = overError + underError;
                }
            }//End of while overflow > 1

            solution = new List<double>();
            foreach (double localFlow in localFlows)
                solution.Add(localFlow);

            underflow = underError;
        }

        private static Dictionary<double, double> BuildVirtualDic(List<double> ratios, List<Dictionary<double, double>> units, int size)
        {
            Dictionary<double, double> virtualMixed = new Dictionary<double, double>(size);
            for (int i = 0; i < ratios.Count; i++)
            {
                foreach (double key in units[i].Keys)
                    if (!virtualMixed.ContainsKey(key))
                        virtualMixed.Add(key, units[i][key] * ratios[i]);
                    else
                        virtualMixed[key] += units[i][key] * ratios[i];
            }
            return virtualMixed;
        }

        private static double FindLocalMaxima(Dictionary<double, double> unit, Dictionary<double, double> mixed)
        {
            double minUnitRatio = double.MaxValue;
            foreach (double key in unit.Keys)
            {
                if (unit[key] > 0)
                {
                    double local = mixed[key] / unit[key];
                    if (local < minUnitRatio)
                        minUnitRatio = local;
                }
            }
            if (minUnitRatio == double.MaxValue)
                return 0;
            else
                return minUnitRatio;
        }
        
        private static double FindLocalMinima(Dictionary<double, double> unit, Dictionary<double, double> mixed, Dictionary<double, double> virtualMixed)
        {
            double minUnitRatio = double.MaxValue;
            foreach (double key in unit.Keys)
            {
                if (unit[key] > 0 && virtualMixed[key] > mixed[key])
                {
                    double local = (virtualMixed[key] - mixed[key]) / unit[key];
                    if (local < minUnitRatio)
                        minUnitRatio = local;
                }
            }
            if (minUnitRatio == double.MaxValue)
                return 0;
            else
                return minUnitRatio;
        }
    }
}
