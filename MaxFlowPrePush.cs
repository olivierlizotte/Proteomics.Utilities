using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteomics.Utilities
{
    public class MaxFlowPrePush
    {
        private static double ComputeOverflow(List<List<double>> fragRatios, List<long> ratios, List<long> fragments)
        {
            double error = 0;
            for (int i = 0; i < fragments.Count; i++)
            {
                double sum = 0;
                for (int j = 0; j < ratios.Count; j++)
                    sum += fragRatios[j][i] * ratios[j];
                if (sum > fragments[i])
                    error += sum - fragments[i];
            }
            return error;
        }

        private static double ComputeUnderflow(List<List<double>> fragRatios, List<long> ratios,
                                        List<long> fragments)
        {
            double error = 0;
            for (int i = 0; i < fragments.Count; i++)
            {
                double sum = 0;
                for (int j = 0; j < ratios.Count; j++)
                    sum += fragRatios[j][i] * ratios[j];
                if (sum < fragments[i])
                    error += fragments[i] - sum;
            }
            return error;
        }

        private static long FindLocalMaximumFlow(List<double> fragRatio, List<long> fragments, long sumOfIntensities)
        {
            long cumul = 1;
            //Check if there is at least one common fragment
            bool isPresent = false;
            for (int i = 0; i < fragments.Count; i++)
                if (fragRatio[i] > 0 && fragments[i] > 0)
                    isPresent = true;
            if (isPresent)
            {
                while (cumul <= sumOfIntensities)
                {
                    for (int i = 0; i < fragments.Count; i++)
                    {
                        if (fragRatio[i] * cumul > fragments[i])
                            return cumul;
                    }
                    cumul++;
                }
            }
            return cumul;
        }

        private static void NormalizeList(List<double> list, double dividend)
        {
            double sum = 0.0;
            foreach (double item in list)
                sum += item;
            sum /= dividend;
            for (int i = 0; i < list.Count; i++)
                list[i] /= sum;
        }

        private static double ComputeMaxFlow(List<List<double>> fragRatios,
                                    List<long> fragments,
                                ref List<List<double>> optimalSolutions)
        {
            //Lists must have same number of fragments, ordered in the same manner
            long sumOfIntensities = 0;
            for (int i = 0; i < fragments.Count; i++)
                sumOfIntensities += fragments[i];

            List<long> localFlows = new List<long>();
            foreach (List<double> fragmentRatio in fragRatios)
                localFlows.Add(FindLocalMaximumFlow(fragmentRatio, fragments, sumOfIntensities));

            double overError = ComputeOverflow(fragRatios, localFlows, fragments);
            double underError = ComputeUnderflow(fragRatios, localFlows, fragments);
            int bestIndex = 0;

            int iterSize = 1;
            double bestOverallError = double.MaxValue;
            List<long> bestLocalFlows = new List<long>();

            while (overError > 0 && iterSize < 10000)
            {
                bestIndex = -1;
                double smallestUnderError = double.MaxValue;
                double smallestOverError = overError;
                for (int i = 0; i < fragRatios.Count; i++)
                {
                    if (localFlows[i] > 0)
                    {
                        localFlows[i] -= iterSize;
                        double tmpErrorMinus = ComputeUnderflow(fragRatios, localFlows, fragments);
                        double tmpErrorPlus  =  ComputeOverflow(fragRatios, localFlows, fragments);

                        if (tmpErrorPlus < overError && (tmpErrorMinus < smallestUnderError 
                            || (tmpErrorMinus == smallestUnderError && tmpErrorPlus < smallestOverError)))
                        {
                            smallestOverError = tmpErrorPlus;
                            smallestUnderError = tmpErrorMinus;
                            bestIndex = i;
                        }
                        localFlows[i] += iterSize;
                    }
                }
                if (bestIndex != -1)
                {
                    localFlows[bestIndex] -= iterSize;
                    iterSize = 1;
                }
                else
                    iterSize++;

                overError = ComputeOverflow(fragRatios, localFlows, fragments);
                underError = ComputeUnderflow(fragRatios, localFlows, fragments);
                if (overError + underError < bestOverallError)
                {
                    bestLocalFlows = new List<long>(localFlows);
                    bestOverallError = overError + underError;
                }
            }
            optimalSolutions.Clear();

            List<double> newList = new List<double>();
            foreach (int localFlow in localFlows)
                newList.Add(localFlow);
            optimalSolutions.Add(newList);

            newList = new List<double>();
            foreach (int localFlow in bestLocalFlows)
                newList.Add(localFlow);
            optimalSolutions.Add(newList);

            return ComputeUnderflow(fragRatios, localFlows, fragments);
        }//*/

        public static List<double> Compute(List<List<double>> ratiosToFit, List<string> ratioNames, int precision, List<double> capacity, ref double overFlow, ref double underFlow, ref double errorInPercent)
        {
            //double 
            /*
            foreach (List<double> list in ratiosToFit)
                NormalizeList(list, precision);//*/

            List<List<double>> solutions = new List<List<double>>();
            List<long> expandedCapacity = new List<long>();
            foreach (double val in capacity)
                expandedCapacity.Add((long)(val*precision));

            double error = ComputeMaxFlow(ratiosToFit, expandedCapacity, ref solutions);

            //Compute average
            List<double> average = new List<double>();
            for (int i = 0; i < solutions[0].Count; i++)
            {
                double sum = 0.0;
                foreach (List<double> solution in solutions)
                    sum += solution[i];
                double avg = sum / (double)solutions.Count;
                average.Add(avg);
            }
            solutions.Add(average);

            //Compute expected error in percentage
            double errorCumul = 0.0;
            for (int k = 0; k < ratiosToFit[0].Count; k++)
            {
                double peakSum = 0;
                for (int i = 0; i < ratiosToFit.Count; i++)
                    peakSum += ratiosToFit[i][k] * average[i];
                errorCumul += Math.Abs(expandedCapacity[k] - peakSum);
            }

            double sumOfIntensities = 0;
            foreach (long item in expandedCapacity)
                sumOfIntensities += item;
            /*
            Console.WriteLine(" -=+ Error cumulated : " + (float)(errorCumul / (double)sumOfIntensities) * 100 + "% +=- ");
            
            foreach (List<double> solution in solutions)
            {
                Console.WriteLine(" ----------------- --------------- ----------------- ------------------ ---------------- ");
                Console.WriteLine("Max Flow computed (error of " + error / (double)(precision * precision) + ")");
                for (int i = 0; i < ratioNames.Count; i++)
                    Console.WriteLine("     " + ratioNames[i] + " -> " + solution[i] / (double) precision + "");
            }
            Console.WriteLine(" ----------------- --------------- ----------------- ------------------ ---------------- ");
            Console.WriteLine("Number of solutions : " + solutions.Count);
            //*/
            List<double> result = new List<double>();
            double sumResultRatios = 0;
            foreach(double val in average)
            {
                sumResultRatios += val / (double)precision;
                result.Add(val / (double)precision);
            }
            Console.WriteLine("Sum of result ratios for average solution : " + sumResultRatios);
            errorInPercent = (errorCumul / sumOfIntensities) * 100.0;
            overFlow = 0;
            underFlow = error;
            return result;
        }
    }
}
