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
    public interface ITargetDecoy
    {
        bool Target { get; }
        bool Decoy { get; }
    }
    ///TODO Wrap the comparison class in a second one, forcing decoys to appear before targets in the case the 
    ///comparison returns '0'

    public class CompareDecoyFirst<T> : IComparer<T> where T : ITargetDecoy
    {
        private Comparison<T> fct;
        public CompareDecoyFirst(Comparison<T> compareFct)
        {
            fct = compareFct;
        }

        // Calls CaseInsensitiveComparer.Compare with the parameters reversed. 
        public int Compare(T x, T y)
        {
            int rez = fct(x, y);
            if (rez == 0)
                return x.Target.CompareTo(y.Target);//-x.Decoy.CompareTo(y.Decoy);
            else
                return rez;
        }
    }

    /// <summary>
    /// This class wraps a List of ITargetDecoy objects. Based on a list of Comparison functions, this tool
    /// reports the biggest list of target it can find at a specified FDR. 
    /// Watch out for overfitting!
    /// TODO : this is only for testing!! This class is used to test the impact of different classifications.
    /// The final FDRizer code will only use one method, optimized to report the most target, without overfiting 
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public class FDRizer<T> where T : ITargetDecoy
    {
        public IterHelper<T> helper;
        public FDRizer(List<T> items, List<Comparison<T>> comparers, T nullObject)
        {
            helper = new IterHelper<T>(items, nullObject);
            foreach (Comparison<T> fct in comparers)
                helper.AddComparison(new CompareDecoyFirst<T>(fct));

            List<T> mashupList = BuildMashUpList(helper.Lists);
            helper.Lists.Add(mashupList);
            helper.Indexes.Add(0);
            helper.ReverseIndexes.Add(mashupList.Count - 1);            
        }

        public class IterHelper<_T> where _T : ITargetDecoy
        {
            public List<int> Indexes;
            public List<int> ReverseIndexes;
            public List<List<_T>> Lists;
            public List<_T> unSortedList;
                        
            Dictionary<_T, bool> cumulTotal;
            List<_T> cumulTarget;
            List<_T> cumulDecoy;
            Dictionary<_T, int> ignoredDecoys;
            /*
            private void SortLists(double desired_fdr)
            {
                List<int> targetsPerList = new List<int>();
                for(int i = 0; i < Lists.Count; i++)
                    targetsPerList.Add(ComputeAtFDR(Lists[i], desired_fdr).Count);

                for(int i = 0; i < Lists.Count; i++)
                    for(int j = i + 1; j < Lists.Count; j++)
                        if (targetsPerList[j] < targetsPerList[i])
                        {
                            List<_T> swap = Lists[i];
                            Lists[i] = Lists[j];
                            Lists[j] = swap;

                            int iSwap = targetsPerList[i];
                            targetsPerList[i] = targetsPerList[j];
                            targetsPerList[j] = iSwap;
                        }
            }//*/

            private _T emptyOrNullObject;

            public IterHelper(List<_T> items, _T nullObject)
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

            public List<_T> TryToDropFDRbkp(Dictionary<_T, int> cumul)
            {
                List<_T> decoys = new List<_T>();
                for (int i = 0; i < Lists.Count; i++)
                    while (ReverseIndexes[i] - 1 >= Indexes[i] && (!Lists[i][ReverseIndexes[i] - 1].Target || !cumul.ContainsKey(Lists[i][ReverseIndexes[i] - 1])))
                    {
                        ReverseIndexes[i]--;
                        decoys.Add(Lists[i][ReverseIndexes[i]]);
                    }
                return decoys;
            }

            public List<_T> GetWorstTargetFromList(Dictionary<_T, int> cumul)
            {
                List<_T> decoys = new List<_T>();
                for (int i = 0; i < Lists.Count; i++)
                    while (ReverseIndexes[i] - 1 >= Indexes[i] && (!Lists[i][ReverseIndexes[i] - 1].Target || !cumul.ContainsKey(Lists[i][ReverseIndexes[i] - 1])))
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
                                                /*
            private List<T> GetListUntilNbDecoy(int indexList, Dictionary<T, int> cumul, int nbCumulDecoy, int nbCumulTarget, int nbDecoysToAdd, ref int bestIndex, ref int addedDecoy)
            {
                int cumulDecoy = 0;
                int cumulTarget = 0;
                bestIndex = Indexes[indexList];
                
                List<T> results = new List<T>();
                List<T> missingResults = new List<T>();
                for (int index = Indexes[indexList]; index < Lists[indexList].Count; index++)
                {
                    T match = Lists[indexList][index];
                    if (!cumul.ContainsKey(match))
                    {
                        if (match.Target)
                        {
                            cumulTarget++;
                            if (missingResults.Count > 0)
                            {
                                results.AddRange(missingResults);
                                missingResults.Clear();
                            }
                            bestIndex = index;
                            results.Add(match);
                        }
                        else
                        {
                            if(cumulDecoy >= nbDecoysToAdd)
                                break;

                            cumulDecoy++;
                            missingResults.Add(match);
                        }
                    }
                }
                addedDecoy = cumulDecoy - missingResults.Count;
                return results;
            }//*/

            private void MatchLastDecoys()
            {
                //Ignore last Decoys
                foreach (_T ignorable in GetLastDecoys())
                {
                    ignoredDecoys.Add(ignorable, 0);
                    if (!cumulTotal.ContainsKey(ignorable))
                        cumulTotal.Add(ignorable, ignorable.Target);
                    else
                        if (ignorable.Target)
                            cumulTarget.Remove(ignorable);
                        else
                            cumulDecoy.Remove(ignorable);
                }
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

            public List<_T> ComputeDifferently_bkp(double desired_fdr)
            {
                Dictionary<_T, bool> result = new Dictionary<_T, bool>();
                foreach (_T match in unSortedList)
                    result.Add(match, match.Target);

                //Compute list of items based on end of lists
                int nbElemAtATime = 1;
                while (ComputeFDR(result) > desired_fdr && nbElemAtATime < unSortedList.Count)
                {
                    int bestNbDecoy = 0;
                    int bestDecoyList = 0;
                    for (int i = 0; i < Lists.Count; i++)
                    {
                        int nbDecoy = GetNbDecoy(i, nbElemAtATime);
                        if (nbDecoy > bestNbDecoy)
                        {
                            bestNbDecoy = nbDecoy;
                            bestDecoyList = i;
                        }
                    }
                    if (bestNbDecoy >= nbElemAtATime)
                    {
                        for (int i = 0; i < nbElemAtATime + bestNbDecoy - 1 && ReverseIndexes[bestDecoyList] > 0; i++)
                        {
                            _T match = Lists[bestDecoyList][ReverseIndexes[bestDecoyList]];
                            if (result.ContainsKey(match))
                                result.Remove(match);
                            ReverseIndexes[bestDecoyList]--;
                        }
                    }
                    else
                        nbElemAtATime++;
                }

                bool keepGoing = true;
                int nbTotalDecoy = 0;
                while (keepGoing)
                {
                    nbTotalDecoy = 0;
                    foreach (_T match in result.Keys)
                        if (match.Decoy)
                            nbTotalDecoy++;
                    //Add stuff to meet desired fdr
                    int bestList = 0;
                    List<_T> bestListResult = new List<_T>();
                    int bestIndex = 0;
                    int realNbDecoyAdded = 0;//*/
                    for (int i = 0; i < Lists.Count; i++)
                    {
                        int addedDecoy = 0;
                        int addedIndex = 0;
                        //List<T> tmpResult = GetFDRedList(i, cumul, nbCumulDecoy, nbCumulTarget, nbToleratedDecoy / (double)this.unSortedList.Count, ref addedIndex, ref addedDecoy);
                        List<_T> tmpResult = GetFDRedList(i, result, nbTotalDecoy, result.Count - nbTotalDecoy, desired_fdr, ref addedIndex, ref addedDecoy);
                        if (tmpResult.Count - addedDecoy > bestListResult.Count - realNbDecoyAdded)
                        {
                            realNbDecoyAdded = addedDecoy;
                            bestListResult = tmpResult;
                            bestList = i;
                            bestIndex = addedIndex;
                        }
                    }

                    if (bestListResult.Count > 0)
                    {
                        //TODO see if it improves results at 1% FDR
                        for (int i = 0; i < (bestListResult.Count + 1) / 2; i++)//bestListResult.Count; i++)//
                        {
                            Indexes[bestList]++;
                            _T match = bestListResult[i];
                            if (!result.ContainsKey(match))
                                result.Add(match, match.Target);
                        }
                    }
                    else
                        keepGoing = false;
                }
                if (result.Count == nbTotalDecoy || nbTotalDecoy / (double)(result.Count - nbTotalDecoy) > desired_fdr)
                    return new List<_T>();
                else
                    return new List<_T>(result.Keys);
            }

            public class iWork { public bool added = false; public bool removed = false; }

            private bool AllDone(Dictionary<_T, iWork> result)
            {
                foreach (iWork val in result.Values)
                    if (!(val.added || val.removed))
                        return false;
                return true;
            }
            public List<_T> ComputeDifferently(double desired_fdr)
            {
                Dictionary<_T, iWork> result = new Dictionary<_T, iWork>();
                foreach(_T match in unSortedList)
                    result.Add(match, new iWork());

                bool keepGoing = true;

                int nbTotalDecoy = 0;
                int nbTotalTarget = 0;
                do
                {                    
                    nbTotalDecoy = 0;
                    nbTotalTarget = 0;
                    foreach (_T match in result.Keys)
                        if (result[match].added)
                            if(match.Decoy)
                                nbTotalDecoy++;
                            else
                                nbTotalTarget++;

                    //Add stuff to meet desired fdr
                    int bestList = 0;                    
                    int bestIndex = 0;
                    int bestReverseIndex = 0;
                    int bestAdvance = 0;//*/
                    for (int k = 0; k < Lists.Count; k++)
                    {
                        int newIndex = 0;
                        int newReverse = 0;
                        
                        int advance = GetIndexes(k, result, desired_fdr, ref newIndex, ref newReverse, nbTotalTarget, nbTotalDecoy);
                        for (int i = 0; i < ((newIndex - Indexes[k]) + 1) / 2; i++)
                        {
                            result[Lists[k][Indexes[k]]].added = true;
                            Indexes[k]++;
                        }
                        for (int i = 0; i < ((ReverseIndexes[k] - newReverse) + 1) / 2; i++)
                        {
                            result[Lists[k][ReverseIndexes[k]]].removed = true;
                            ReverseIndexes[k]--;
                        }//*/

                        if (advance > bestAdvance)
                        {
                            bestAdvance = advance;
                            bestList = k;
                            bestIndex = newIndex;
                            bestReverseIndex = newReverse;
                        }
                    }/*
                    //if(bestIndex > Indexes[bestList])//bestAdvance > 0)
                    //{
                        //if ((bestIndex - Indexes[bestList]) > (ReverseIndexes[bestList] - bestReverseIndex))
                        //{
                            int maxI = ((bestIndex - Indexes[bestList]) + 1) / 2;
                            for (int i = 0; i < maxI; i++)
                            {
                                result[Lists[bestList][Indexes[bestList]]].added = true;
                                Indexes[bestList]++;
                            }
                    //}
                    //if(bestReverseIndex < ReverseIndexes[bestList])
                    //{
                        //else
                        //{//Beat 7401
                            maxI = ((ReverseIndexes[bestList] - bestReverseIndex) + 1) / 2;
                            for (int i = 0; i < maxI; i++)
                            {
                                result[Lists[bestList][ReverseIndexes[bestList]]].removed = true;
                                ReverseIndexes[bestList]--;
                            }
                        //}//*/
                    //}
                    if(bestAdvance == 0)
                        keepGoing = false;//*/
                }
                while(keepGoing);//!AllDone(result));                

                List<_T> output = new List<_T>();
                foreach(_T key in result.Keys)
                    if(result[key].added && !result[key].removed)
                        output.Add(key);

                if (nbTotalTarget > 0 && nbTotalDecoy / (double) nbTotalTarget <= desired_fdr)
                    return output;
                else
                    return new List<_T>();
                
                /*

                Dictionary<T, bool> result = new Dictionary<T, bool>();                
                foreach (T match in unSortedList)
                    result.Add(match, match.Target);
                                
                while (ComputeFDR(result) > desired_fdr)
                {
                    int nbTotalDecoy = 0;
                    foreach (T match in result.Keys)
                        if (match.Decoy)
                            nbTotalDecoy++;
                    int nbElemAtATime = (int) (1 + (result.Count - nbTotalDecoy) * desired_fdr);
                    int bestNbDecoy = 0;
                    int bestList = 0;
                    for (int i = 0; i < Lists.Count; i++)
                    {
                        int nbDecoy = GetNbDecoy(i, nbElemAtATime);
                        if (nbDecoy > bestNbDecoy)
                        {
                            bestNbDecoy = nbDecoy;
                            bestList = i;
                        }
                    }
                    if (bestNbDecoy > 0)
                    {
                        for (int i = 0; i < (nbElemAtATime + 1) / 2; i++)
                        {
                            T match = Lists[bestList][ReverseIndexes[bestList]];
                            if (result.ContainsKey(match))
                                result.Remove(match);
                            ReverseIndexes[bestList]--;
                        }
                    }
                    else
                    {//*/
                        /*
                        int bestNbTarget = 0;
                        for (int i = 0; i < Lists.Count; i++)
                        {
                            int nbTarget = GetNbTarget(i, nbElemAtATime);
                            if (nbTarget > bestNbTarget)
                            {
                                bestNbTarget = nbTarget;
                                bestList = i;
                            }
                        }
                        if (bestNbTarget > 0)
                        {
                            for (int i = 0; i < nbElemAtATime && Indexes[bestList] <= ReverseIndexes[bestList]; i++)
                            {
                                T match = Lists[bestList][Indexes[bestList]];
                                if (!result.ContainsKey(match))
                                    result.Add(match, match.Target);
                                Indexes[bestList]++;
                            }
                        }
                        else//*/
                     //   nbElemAtATime++;
                   // }
                //}
                //return new List<T>(result.Keys);
            }

            public List<_T> ComputeRtoL(double desired_fdr)
            {
                Dictionary<_T, bool> result = new Dictionary<_T, bool>();
                foreach (_T match in unSortedList)
                    result.Add(match, match.Target);

                int nbElemAtATime = 1;
                while (ComputeFDR(result) > desired_fdr && nbElemAtATime < unSortedList.Count)
                {
                    int bestNbDecoy = 0;
                    int bestList = 0;
                    for (int i = 0; i < Lists.Count; i++)
                    {
                        int nbDecoy = GetNbDecoy(i, nbElemAtATime);
                        if (nbDecoy > bestNbDecoy)
                        {
                            bestNbDecoy = nbDecoy;
                            bestList = i;
                        }
                    }
                    if (bestNbDecoy > 0)
                    {
                        for (int i = 0; i < nbElemAtATime + bestNbDecoy - 1 && ReverseIndexes[bestList] >= 0; i++)
                        {
                            _T match = Lists[bestList][ReverseIndexes[bestList]];
                            if (result.ContainsKey(match))
                                result.Remove(match);
                            ReverseIndexes[bestList]--;
                        }
                    }
                    else
                    {/*
                        int bestNbTarget = 0;
                        for (int i = 0; i < Lists.Count; i++)
                        {
                            int nbTarget = GetNbTarget(i, nbElemAtATime);
                            if (nbTarget > bestNbTarget)
                            {
                                bestNbTarget = nbTarget;
                                bestList = i;
                            }
                        }
                        if (bestNbTarget > 0)
                        {
                            for (int i = 0; i < nbElemAtATime && Indexes[bestList] <= ReverseIndexes[bestList]; i++)
                            {
                                T match = Lists[bestList][Indexes[bestList]];
                                if (!result.ContainsKey(match))
                                    result.Add(match, match.Target);
                                Indexes[bestList]++;
                            }
                        }
                        else//*/
                            nbElemAtATime++;
                    }
                }
                if (ComputeFDR(result) <= desired_fdr)
                    return new List<_T>(result.Keys);
                else
                    return new List<_T>();
            }

            public List<_T> ComputeLtoR(double desired_fdr)
            {
                //Step 1: Sort lists in descending order of impact (number of targets at given fdr)
                //SortLists(desired_fdr);
                
                double reverseFactor = 1.0;

                do
                {
                    //Ignore last Decoys
                    //MatchLastDecoys();

                    int worstList = 0;
                    List<_T> worstListResult = new List<_T>();
                    int worstIndex = 0;
                    int realNbTargetAdded = 0;
                    for (int i = 0; i < Lists.Count; i++)
                    {
                        int addedDecoy = 0;
                        int addedIndex = 0;
                        List<_T> tmpResult = GetReverseFDRedList(i, cumulTotal, ignoredDecoys, cumulDecoy.Count, cumulTarget.Count, desired_fdr * reverseFactor, ref addedIndex, ref addedDecoy);
                        if (addedDecoy > realNbTargetAdded)
                        {
                            realNbTargetAdded = addedDecoy;
                            worstListResult = tmpResult;
                            worstList = i;
                            worstIndex = addedIndex;
                        }
                    }
                    if (worstListResult.Count > 0)
                    {
                        ReverseIndexes[worstList] = worstIndex;
                        //worstListResult.AddRange(GetLastDecoys());
                        foreach (_T ignorable in worstListResult)
                        {
                            if (!ignoredDecoys.ContainsKey(ignorable))
                                ignoredDecoys.Add(ignorable, 0);

                            if (!cumulTotal.ContainsKey(ignorable))
                                cumulTotal.Add(ignorable, ignorable.Target);
                            else
                                if (ignorable.Target)
                                    cumulTarget.Remove(ignorable);
                                else
                                    cumulDecoy.Remove(ignorable);
                        }
                    }
                    else
                    {
                        int bestList = 0;
                        List<_T> bestListResult = new List<_T>();
                        int bestIndex = 0;
                        int realNbDecoyAdded = 0;//*/
                        for (int i = 0; i < Lists.Count; i++)
                        {
                            int addedDecoy = 0;
                            int addedIndex = 0;
                            //List<T> tmpResult = GetFDRedList(i, cumul, nbCumulDecoy, nbCumulTarget, nbToleratedDecoy / (double)this.unSortedList.Count, ref addedIndex, ref addedDecoy);
                            List<_T> tmpResult = GetFDRedList(i, cumulTotal, cumulDecoy.Count, cumulTarget.Count, desired_fdr, ref addedIndex, ref addedDecoy);
                            if (tmpResult.Count - addedDecoy > bestListResult.Count - realNbDecoyAdded)
                            {
                                realNbDecoyAdded = addedDecoy;
                                bestListResult = tmpResult;
                                bestList = i;
                                bestIndex = addedIndex;
                            }
                        }

                        if (bestListResult.Count > 0)
                        {
                            //TODO see if it improves results at 1% FDR
                            for (int i = 0; i < (bestListResult.Count + 1) / 2; i++)//bestListResult.Count; i++)//
                            {
                                Indexes[bestList]++;
                                _T precursor = bestListResult[i];
                                if (!cumulTotal.ContainsKey(precursor))
                                {
                                    if (precursor.Target)
                                        cumulTarget.Add(precursor);
                                    else
                                        cumulDecoy.Add(precursor);

                                    cumulTotal.Add(precursor, precursor.Target);
                                }
                            }
                        }
                        else
                        {/*
                            List<T> toAdd = GetNextTargets(cumulTotal);
                            if (toAdd.Count > 0)
                            {
                                foreach (T precursor in toAdd)
                                {
                                    if (!cumulTotal.ContainsKey(precursor))
                                    {
                                        if (precursor.Target)
                                            cumulTarget.Add(precursor);
                                        else
                                            cumulDecoy.Add(precursor);

                                        cumulTotal.Add(precursor, bestList);
                                    }
                                }
                            }
                            else//*/
                                reverseFactor++;
                        }
                    }
                } while (reverseFactor * desired_fdr < 0.25 && reverseFactor < 100);
                foreach (_T ignorable in ignoredDecoys.Keys)
                    cumulTotal.Remove(ignorable);

                if (ComputeFDR(cumulTotal) <= desired_fdr)
                    return new List<_T>(cumulTotal.Keys);
                else
                    return new List<_T>();
            }

            private List<_T> GetReverseFDRedList(int indexList, Dictionary<_T, bool> cumul, Dictionary<_T, int> ignoredDecoys, int nbCumulDecoy, int nbCumulTarget, double desired_fdr, ref int bestIndex, ref int addedDecoy)
            {
                int cumulDecoy = 0;
                int cumulTarget = 0;
                bestIndex = Indexes[indexList];

                List<_T> results = new List<_T>();
                List<_T> missingResults = new List<_T>();
                for (int index = ReverseIndexes[indexList] - 1; index >= Indexes[indexList]; index--)
                {
                    _T match = Lists[indexList][index];
                    //if (!cumul.ContainsKey(match))
                    {
                        if (match.Target)
                        {
                            if (cumul.ContainsKey(match) && !ignoredDecoys.ContainsKey(match))
                                break;
                            else
                                cumulTarget++;
                        }
                        else
                            cumulDecoy++;
                        if (!match.Target && cumulTarget / (double)cumulDecoy <= desired_fdr*20)//TODO Check if its still efficient
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

            private int GetIndexes(int indexList, Dictionary<_T, iWork> cumul, double desired_fdr, ref int bestIndex, ref int bestReverseIndex, int cumulTarget, int cumulDecoy)
            {
                int iRev = ReverseIndexes[indexList];
                int iInd = Indexes[indexList];
                int addedDecoy = 0;
                int removedDecoy = 0;
                int addedTarget = 0;
                int removedTarget = 0;
                int bestNbAddedTarget = 0;
                bestReverseIndex = ReverseIndexes[indexList];
                bestIndex = Indexes[indexList];
                while (iRev > iInd)
                {
                    if (addedTarget - removedTarget >= bestNbAddedTarget && (cumulDecoy + addedDecoy - removedDecoy) / (double)(cumulTarget + addedTarget - removedTarget) <= desired_fdr)
                    {
                        bestNbAddedTarget = addedTarget - removedTarget;
                        bestIndex = iInd;
                        bestReverseIndex = iRev;
                    }

                    _T matchInd = Lists[indexList][iInd];
                    if (matchInd.Target)
                    {
                        if (!cumul[matchInd].added && !cumul[matchInd].removed)
                            addedTarget++;
                        iInd++;
                    }
                    else
                    {
                        _T matchRev = Lists[indexList][iRev];
                        if (matchRev.Decoy)
                        {
                            if (!cumul[matchRev].removed && cumul[matchRev].added)
                                removedDecoy++;
                            iRev--;
                        }
                        else
                        {
                            if (!cumul[matchInd].added && !cumul[matchInd].removed)
                                addedDecoy++;
                            if (cumul[matchRev].added && !cumul[matchRev].removed)
                                removedTarget++;
                            iInd++;
                            iRev--;
                        }
                    }
                }
                return ReverseIndexes[indexList] - bestReverseIndex + bestIndex - Indexes[indexList];
            }
        }
        
        public class iWork
        {
            public bool added = false;
            public bool removed = false;
            public T match;
            public double score = 0;
            public iWork(T pMatch)
            {
                this.match = pMatch;
            }
            public static int CompareDescendingScore(iWork left, iWork right)
            {
                int compare = -left.score.CompareTo(right.score);
                if (compare == 0)
                    return left.match.Target.CompareTo(right.match.Target);//Decoy first, target after
                else
                    return compare;
            }
        }

        private List<T> ComputeFDRedList(List<iWork> orderedList, double desired_fdr)
        {
            int cumulDecoy = 0;
            int cumulTarget = 0;

            List<T> results = new List<T>();
            List<T> missingResults = new List<T>();
            for (int index = 0; index < orderedList.Count; index++)
            {
                iWork item = orderedList[index];

                if (item.match.Target)
                    cumulTarget++;
                else
                    cumulDecoy++;

                if (item.match.Target && cumulDecoy / (double)cumulTarget <= desired_fdr)
                {
                    if (missingResults.Count > 0)
                    {
                        results.AddRange(missingResults);
                        missingResults.Clear();
                    }
                    results.Add(item.match);
                }
                else
                    missingResults.Add(item.match);
            }
            return results;
        }
        
        public List<T> BuildMashUpList(List<List<T>> Lists)
        {
            //Compute Scores
            Dictionary<T, iWork> scores = new Dictionary<T, iWork>();
            foreach (T match in Lists[0])
                if (!scores.ContainsKey(match)) 
                    scores.Add(match, new iWork(match));

            for (int indexList = 0; indexList < Lists.Count; indexList++)
            {
                int cumulDecoy = 0;
                for (int i = 0; i < Lists[indexList].Count; i++)
                    if (scores[Lists[indexList][i]].match.Decoy)
                        cumulDecoy++;
                double totalDecoy = cumulDecoy;

                cumulDecoy = 0;
                for (int i = 0; i < Lists[indexList].Count; i++)
                {
                    iWork item = scores[Lists[indexList][i]];
                    if (item.match.Decoy)
                        cumulDecoy++;
                    item.score += 1.0 - cumulDecoy / totalDecoy;
                }
            }

            //Sort based on score
            List<iWork> listScores = new List<iWork>(scores.Values);
            listScores.Sort(iWork.CompareDescendingScore);

            List<T> sortedList = new List<T>();
            foreach (iWork item in listScores)
                sortedList.Add(item.match);
            return sortedList;             
        }

        public static List<T> ComputeAtFDR(List<T> list, double desired_fdr)
        {
            int cumulDecoy = 0;
            int cumulTarget = 0;

            List<T> results = new List<T>();
            List<T> missingResults = new List<T>();
            for (int index = 0; index < list.Count; index++)
            {
                if (list[index].Target)
                    cumulTarget++;
                else
                    cumulDecoy++;

                if (list[index].Target && cumulDecoy / (double)cumulTarget <= desired_fdr)
                {
                    if (missingResults.Count > 0)
                    {
                        results.AddRange(missingResults);
                        missingResults = new List<T>();
                    }
                    results.Add(list[index]);
                }
                else
                    missingResults.Add(list[index]);
            }
            return results;
        }//*/

        public List<T> FDRByProbability(double desiredFDR)
        {
            List<T> probList = BuildMashUpList(helper.Lists);
            return ComputeAtFDR(probList, desiredFDR);
        }

        public void ReStart()
        {
            helper.ReStart();
        }

        public List<T> Launch(double desired_fdr, bool displayValues = false)
        {
            //if(displayValues)
            //    Sol.CONSOLE.OutputLine("Optimizing scores to get maximum precursors at " + (100 * desired_fdr) + "% fdr");

            List<T> result = helper.ComputeRtoL(desired_fdr);
            helper.ReStart();
            List<T> resultA = helper.ComputeLtoR(desired_fdr);
            if (resultA.Count > result.Count)
                result = resultA;
            helper.ReStart();
            resultA = FDRByProbability(desired_fdr);
            if (resultA.Count > result.Count)
                result = resultA;
            //resultA = helper.ComputeDifferently_bkp(desired_fdr);
            helper.ReStart();
            resultA = helper.ComputeDifferently(desired_fdr);
            if (resultA.Count > result.Count)
                result = resultA;
            int targets = 0;
            int decoys = 0;
            foreach (T elem in result)
            {
                if (elem.Target)
                    targets++;
                if (elem.Decoy)
                    decoys++;
            }
            if (targets <= 0 || decoys / (double)targets > desired_fdr)
                result = new List<T>();

            if (displayValues)
            {
                Console.WriteLine(" === Score impact distribution === ");
                for (int i = 0; i < helper.Indexes.Count; i++)
                {
                    double ratio = helper.Indexes[i] / (double)helper.unSortedList.Count;
                    double reverseRatio = helper.ReverseIndexes[i] / (double)helper.unSortedList.Count;
                    Console.WriteLine(">    " + i + ": " + ratio + " (" + reverseRatio + ")");
                }//*/
            }
            return result;
        }
    }
}
