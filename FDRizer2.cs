/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteomics.Utilities
{
    /// <summary>
    /// This class wraps a List of ITargetDecoy objects. Based on a list of Comparison functions, this tool
    /// reports the biggest list of target it can find at a specified FDR. 
    /// Watch out for overfitting!
    /// TODO : this is only for testing!! This class is used to test the impact of different classifications.
    /// The final FDRizer code will only use one method, optimized to report the most target, without overfiting 
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public class FDRizer2<T> where T : ITargetDecoy
    {
        public IterHelper2<T> helper;
        public FDRizer2(List<T> items, List<Comparison<T>> comparers, T nullObject)
        {
            helper = new IterHelper2<T>(items, nullObject);
            foreach (Comparison<T> fct in comparers)
                helper.AddComparison(new CompareDecoyFirst<T>(fct));
        }

        public class IterHelper2<_T> where _T : ITargetDecoy
        {
            public List<int> Indexes;
            public List<int> ReverseIndexes;
            public List<List<_T>> Lists;
            public List<_T> unSortedList;
                        
            Dictionary<_T, bool> cumulTotal;
            List<_T> cumulTarget;
            List<_T> cumulDecoy;
            Dictionary<_T, int> ignoredDecoys;

            private _T emptyOrNullObject;

            public IterHelper2(List<_T> items, _T nullObject)
            {
                emptyOrNullObject = nullObject;
                unSortedList = items;
                Indexes = new List<int>();
                ReverseIndexes = new List<int>();
                Lists = new List<List<_T>>();
                
                cumulTotal = new Dictionary<_T, bool>();
                cumulTarget = new List<_T>();
                cumulDecoy  = new List<_T>();
                ignoredDecoys = new Dictionary<_T, int>();
            }

            public void ReStart()
            {
                for (int i = 0; i < Indexes.Count; i++)
                {
                    Indexes[i] = 0;
                    ReverseIndexes[i] = unSortedList.Count - 1;
                }
                cumulTotal.Clear();
                cumulTarget.Clear();
                cumulDecoy.Clear();
                ignoredDecoys.Clear();
            }
            
            public void AddComparison(IComparer<_T> fct)
            {
                List<_T> list = new List<_T>(unSortedList);
                list.Sort(fct);
                Indexes.Add(0);
                ReverseIndexes.Add(list.Count - 1);
                Lists.Add(list);
            }

            public List<_T> GetLastDecoys()
            {
                List<_T> decoys = new List<_T>();
                for (int i = 0; i < Lists.Count; i++)
                    while (ReverseIndexes[i] - 1 >= Indexes[i] && !Lists[i][ReverseIndexes[i] - 1].Target)
                    {
                        ReverseIndexes[i]--;
                        decoys.Add(Lists[i][ReverseIndexes[i]]);
                    }
                return decoys;
            }

            private int GetNbTarget(int indexList, int nbCycle)
            {
                int nbTarget = 0;
                for (int i = Indexes[indexList]; nbCycle > 0 && i < ReverseIndexes[indexList]; i++)
                {
                    if (Lists[indexList][i].Target)
                        nbTarget++;
                    else                        
                        nbCycle--;
                }
                return nbTarget;
            }

            private int GetNbDecoy(int indexList, int nbCycle)
            {
                int nbDecoy = 0;
                for (int i = ReverseIndexes[indexList]; nbCycle > 0 && i > Indexes[indexList]; i--)
                {
                    if (Lists[indexList][i].Decoy)
                        nbDecoy++;
                    else
                        nbCycle--;
                }
                return nbDecoy;
            }

            public List<_T> GetNextTargets(Dictionary<_T, int> cumul)
            {
                List<_T> targets = new List<_T>();
                for(int i = 0; i < Lists.Count; i++)
                    while(Indexes[i] + 1 < ReverseIndexes[i] && (cumul.ContainsKey(Lists[i][Indexes[i] + 1]) || Lists[i][Indexes[i] + 1].Target))
                    {
                        Indexes[i]++;
                        targets.Add(Lists[i][Indexes[i]]);
                    }
                return targets;
            }                 
            
            private static double ComputeFDR(Dictionary<_T, bool> items)
            {
                int cumulDecoy = 0;
                int cumulTarget = 0;
                foreach(bool isTarget in items.Values)
                    if(isTarget)
                        cumulTarget++;
                    else
                        cumulDecoy++;
                if(cumulTarget > 0)
                    return cumulDecoy / (double)cumulTarget;
                else
                    return 0;
            }

            private List<_T> GetFDRedList(int indexList, Dictionary<_T, bool> cumul, int nbCumulDecoy, int nbCumulTarget, double desired_fdr, ref int bestIndex, ref int addedDecoy)
            {
                int cumulDecoy = 0;
                int cumulTarget = 0;
                bestIndex = Indexes[indexList];

                List<_T> results = new List<_T>();
                List<_T> missingResults = new List<_T>();
                for (int index = bestIndex; index < ReverseIndexes[indexList]; index++)
                {
                    _T match = Lists[indexList][index];
                    if (!cumul.ContainsKey(match))
                    {
                        if (match.Target)
                            cumulTarget++;
                        else
                            cumulDecoy++;
                        if (match.Target && (nbCumulDecoy + cumulDecoy) / (double)(nbCumulTarget + cumulTarget) <= desired_fdr)
                        {
                            if (missingResults.Count > 0)
                            {
                                results.AddRange(missingResults);
                                missingResults.Clear();
                            }
                            bestIndex = index;
                            results.Add(match);
                            addedDecoy = cumulDecoy;
                        }
                        else
                            missingResults.Add(match);
                    }
                }
                return results;
            }
        }

        public List<T> ComputeAtFDR(double desired_fdr)
        {
            Dictionary<T, bool> results = new Dictionary<T, bool>();

            List<List<T>> missingResultsArray = new List<List<T>>();
            int[] missingDecoy = new int[helper.unSortedList.Count];
            int[] missingTarget = new int[helper.unSortedList.Count];
            for (int index = 0; index < helper.Lists.Count; index++)
                missingResultsArray.Add(new List<T>());

            int totalTarget = 0;
            int totalDecoy = 0;

            for (int index = 0; index < helper.unSortedList.Count; index++)
            {
                for (int indexList = 0; indexList < helper.Lists.Count; indexList++)
                {
                    List<T> missingResults = missingResultsArray[indexList];
                    List<T> list = helper.Lists[indexList];
                    if (!results.ContainsKey(list[index]))
                    {
                        if (list[index].Target)
                            missingTarget[indexList]++;
                        else
                            missingDecoy[indexList]++;

                        if (list[index].Target && (totalDecoy + missingDecoy[indexList]) / (double)(missingTarget[indexList] + totalTarget) < desired_fdr)
                        {
                            if (missingResults.Count > 0)
                            {
                                foreach (T item in missingResults)
                                    if (!results.ContainsKey(item))
                                    {
                                        results.Add(item, item.Target);
                                        if (item.Target)
                                            totalTarget++;
                                        else
                                            totalDecoy++;
                                    }
                                missingResults.Clear();
                            }

                            if (!results.ContainsKey(list[index]))
                            {
                                results.Add(list[index], list[index].Target);
                                if (list[index].Target)
                                    totalTarget++;
                                else
                                    totalDecoy++;
                            }
                            missingDecoy[indexList] = 0;
                            missingTarget[indexList] = 0;
                        }
                        else
                            missingResults.Add(list[index]);
                    }
                }
            }
            return results.Keys.ToList<T>();
        }//*/
        
        public void ReStart()
        {
            helper.ReStart();
        }

        public List<T> Launch(double desired_fdr)
        {
            //if(displayValues)
            //    Sol.CONSOLE.OutputLine("Optimizing scores to get maximum precursors at " + (100 * desired_fdr) + "% fdr");

            List<T> result = ComputeAtFDR(desired_fdr);
            int targets = 0;
            int decoys = 0;
            foreach (T elem in result)
            {
                if (elem.Target)
                    targets++;
                if (elem.Decoy)
                    decoys++;
            }

            return result;
        }
    }
}
