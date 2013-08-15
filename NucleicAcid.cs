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
    /// Convertion between Nucleic Acids and Amino Acids (in mamalians)
    /// TODO Make calls more convivial
    /// </summary>
    public static class NucleicAcid
    {
        public static Hashtable AAToNA3;
        public static Hashtable NA3ToAA;
        public static Hashtable AAToTric;
        public static Hashtable TricToAA;

        public static void InitHash()
        {
            //  --  AAToNA  --
            AAToNA3 = new Hashtable();
            AAToNA3["A"] = new string[] { "GCU", "GCC", "GCA", "GCG" };
            AAToNA3["R"] = new string[] { "CGU", "CGC", "CGA", "CGG", "AGA", "AGG" };
            AAToNA3["N"] = new string[] { "AAU", "AAC" };
            AAToNA3["D"] = new string[] { "GAU", "GAC" };
            AAToNA3["C"] = new string[] { "UGU", "UGC" };
            AAToNA3["Q"] = new string[] { "CAA", "CAG" };
            AAToNA3["E"] = new string[] { "GAA", "GAG" };
            AAToNA3["G"] = new string[] { "GGU", "GGC", "GGA", "GGG" };
            AAToNA3["H"] = new string[] { "CAU", "CAC" };
            AAToNA3["I"] = new string[] { "AUU", "AUC", "AUA" };
            AAToNA3["L"] = new string[] { "UUA", "UUG", "CUU", "CUC", "CUA", "CUG" };
            AAToNA3["K"] = new string[] { "AAA", "AAG" };
            AAToNA3["M"] = new string[] { "AUG" };
            AAToNA3["F"] = new string[] { "UUU", "UUC" };
            AAToNA3["P"] = new string[] { "CCU", "CCC", "CCA", "CCG" };
            AAToNA3["S"] = new string[] { "UCU", "UCC", "UCA", "UCG", "AGU", "AGC" };
            AAToNA3["T"] = new string[] { "ACU", "ACC", "ACA", "ACG" };
            AAToNA3["W"] = new string[] { "UGG" };
            AAToNA3["Y"] = new string[] { "UAU", "UAC" };
            AAToNA3["V"] = new string[] { "GUU", "GUC", "GUA", "GUG" };
            AAToNA3["X"] = new string[] { "UAA", "UGA", "UAG" };
            List<string> cumul = new List<string>();
            foreach (string key in AAToNA3.Keys)
                foreach (string na in (string[])AAToNA3[key])
                    cumul.Add(na);
            AAToNA3["*"] = cumul.ToArray();

            //  --  NAToAA  --
            NA3ToAA = new Hashtable();
            foreach (string key in AAToNA3.Keys)
                if (string.Compare(key, "*") != 0)
                    foreach (string na in (string[])AAToNA3[key])
                    {
                        NA3ToAA[na] = key;
                        //NAToAA[na.Replace('U', 'T')] = key;
                    }

            //  --  AAToTric  --
            AAToTric = new Hashtable();
            AAToTric["A"] = "ala";
            AAToTric["R"] = "arg";
            AAToTric["N"] = "asn";
            AAToTric["D"] = "asp";
            AAToTric["C"] = "cys";
            AAToTric["Q"] = "gln";
            AAToTric["E"] = "glu";
            AAToTric["G"] = "gly";
            AAToTric["H"] = "his";
            AAToTric["I"] = "ile";
            AAToTric["L"] = "leu";
            AAToTric["K"] = "lys";
            AAToTric["M"] = "met";
            AAToTric["F"] = "phe";
            AAToTric["P"] = "pro";
            AAToTric["S"] = "ser";
            AAToTric["T"] = "thr";
            AAToTric["W"] = "trp";
            AAToTric["Y"] = "tyr";
            AAToTric["V"] = "val";
            AAToTric["U"] = "sec";
            AAToTric["O"] = "pyl";
            AAToTric["B"] = "asx";
            AAToTric["Z"] = "glx";
            AAToTric["J"] = "xle";
            AAToTric["*"] = "xaa";
            AAToTric["X"] = "ter";

            //  --  TricToAA  --
            TricToAA = new Hashtable();
            foreach (string key in AAToTric.Keys)
                TricToAA[AAToTric[key]] = key;
        }

        public static IEnumerable<char> GetNAs(string naSequence, int index)
        {
            switch (naSequence[index])
            {
                case 'A' : 
                    yield return 'A';
                    break;
                case 'B' :
                    yield return 'G';
                    yield return 'U';
                    yield return 'C';
                    break;
                case 'C':
                    yield return 'C';
                    break;
                case 'D':
                    yield return 'G';
                    yield return 'A';
                    yield return 'U';
                    break;
                case 'G':
                    yield return 'G';
                    break;
                case 'H':
                    yield return 'A';
                    yield return 'C';
                    yield return 'U';
                    break;
                case 'K':
                    yield return 'G';
                    yield return 'U';
                    break;
                case 'M':
                    yield return 'A';
                    yield return 'C';
                    break;
                case 'N':
                    yield return 'A';
                    yield return 'G';
                    yield return 'C';
                    yield return 'U';
                    break;
                case 'R':
                    yield return 'G';
                    yield return 'A';
                    break;
                case 'S':
                    yield return 'G';
                    yield return 'C';
                    break;
                case 'T':
                    yield return 'U';
                    break;
                case 'V':
                    yield return 'G';
                    yield return 'C';
                    yield return 'A';
                    break;
                case 'W':
                    yield return 'A';
                    yield return 'U';
                    break;
                case 'Y':
                    yield return 'U';
                    yield return 'C';
                    break;
            }
        }

        public static IEnumerable<string> ConvertNA3ToSingleAA(string naSequence, int index)
        {
            if(index + 3 < naSequence.Length)
            {
                foreach(char na1 in GetNAs(naSequence, index))
                    foreach(char na2 in GetNAs(naSequence, index + 1))
                        foreach (char na3 in GetNAs(naSequence, index + 2))
                        {
                            string naTriple = na1.ToString() + na2.ToString() + na3.ToString();
                            yield return (string)NA3ToAA[naTriple];
                        }
            }
        }

        public static List<string> ConvertNA3ToAAs(string naSequence)
        {
            List<string> sequences = new List<string>();
            sequences.Add("");
            for (int i = 0; i + 3 < naSequence.Length; i += 3)
            {
                List<string> newList = new List<string>();
                foreach (string aa in ConvertNA3ToSingleAA(naSequence, i))
                    foreach (string poss in sequences)
                        newList.Add(poss + aa);
                if(newList.Count > 0)
                    sequences = newList;
            }
            return sequences;
        }

        public static string ConvertNA3ToAA(string naSequence)
        {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i + 3 < naSequence.Length; i += 3)
            {
                string str = (string)NA3ToAA[naSequence.Substring(i, 3)];
                if (string.IsNullOrEmpty(str))
                    Console.WriteLine("FFFF");
                else
                    sb.Append(str);
            }
            return sb.ToString();
        }

        public static List<string> ConvertAAToNA3s(string aaSequence)
        {
            List<string> possibilities = new List<string>((string[])AAToNA3[aaSequence[0].ToString()]);
            for (int i = 1; i < aaSequence.Length; i++)
            {
                List<string> newList = new List<string>();
                foreach (string item in (string[])AAToNA3[aaSequence[i].ToString()])
                    foreach (string cumulItem in possibilities)
                        newList.Add(cumulItem + item);
                possibilities = newList;
            }
            return possibilities;
        }

        public static bool Test_HASH()
        {
            InitHash();
            string[] nas = new string[] { "A", "C", "G", "U" };
            foreach (string a in nas)
                foreach (string b in nas)
                    foreach (string c in nas)
                        if (NA3ToAA[a + b + c] == null)
                            return false;
            return true;
        }
    }
}
