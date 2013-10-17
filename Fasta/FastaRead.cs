﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace Proteomics.Utilities.Fasta
{
    public static class FastaRead
    {
        public static string Reverse(string str)
        {
            char[] sequence_array = str.ToCharArray();
            Array.Reverse(sequence_array);
            return new string(sequence_array);
        }

        public static IEnumerable<string[]> GetSequences(string fastaFileIn)
        {
            FileStream fs;
            try
            {
                fs = new FileStream(fastaFileIn, FileMode.Open, FileAccess.Read, FileShare.Read);
            }
            catch (System.Exception)
            {
                fs = new FileStream(fastaFileIn, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
            }
                
            using (StreamReader sr = new StreamReader(fs))
            {
                string line;
                string header = null;
                string sequence = "";
                while ((line = sr.ReadLine()) != null)
                {
                    if (line.StartsWith(">"))
                    {
                        if(!string.IsNullOrEmpty(header))
                            yield return new string[]{header, sequence};
                        header = line;
                        sequence = "";
                    }
                    else
                        sequence += line;
                }
                if (!string.IsNullOrEmpty(header))
                    yield return new string[] { header, sequence };
            }
            fs.Close();
        }

        public static void FilterFromIDs(string fastaFileIn, string fastaFileOut, string csvFile)
        {
            vsCSV csv = new vsCSV(csvFile);
            try
            {
                FileStream fs;
                try
                {
                    fs = new FileStream(fastaFileIn, FileMode.Open, FileAccess.Read, FileShare.Read);
                }
                catch (System.Exception)
                {
                    fs = new FileStream(fastaFileIn, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
                }

                FileStream fsOut = new FileStream(fastaFileOut, FileMode.CreateNew);
                using (StreamWriter sw = new StreamWriter(fsOut))
                {
                    using (StreamReader sr = new StreamReader(fs))
                    {
                        string line;
                        bool ignore = true;
                        while ((line = sr.ReadLine()) != null)
                        {
                            if (line.StartsWith(">"))
                            {
                                ignore = true;
                                foreach(string id in csv.LINES_LIST)
                                    if (line.Contains(id))
                                    {
                                        ignore = false;
                                        break;
                                    }
                            }
                            if(!ignore)
                                sw.WriteLine(line);
                        }
                    }
                }
                fsOut.Close();
                fs.Close();
            }
            catch (System.Exception ex)
            {
                Console.WriteLine(ex.Message);
                Console.WriteLine(ex.StackTrace);
            }
        }

        public static void AddIPIToUbiPredFile(string fastaFile, string csvUbiFile, string csvFileOut)
        {
            vsCSV csv = new vsCSV(csvUbiFile);
            vsCSVWriter writer = new vsCSVWriter(csvFileOut);
            try
            {
                FileStream fs;
                try
                {
                    fs = new FileStream(fastaFile, FileMode.Open, FileAccess.Read, FileShare.Read);
                }
                catch (System.Exception)
                {
                    fs = new FileStream(fastaFile, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
                }

                Dictionary<int, string> IDs = new Dictionary<int, string>();
                using (StreamReader sr = new StreamReader(fs))
                {
                    int idNb = 0;
                    string line;
                    while ((line = sr.ReadLine()) != null)
                    {
                        if (line.StartsWith(">"))
                        {
                            idNb++;
                            foreach (string csvLine in csv.LINES_LIST)
                            {
                                string[] splits = csvLine.Split(vsCSV._Generic_Separator);
                                if (line.Contains(splits[1]))
                                    writer.AddLine(csvLine + "," + line.Substring(splits[1].Length + 1, 13));
                            }
                        }
                    }
                }
                fs.Close();
                writer.writeToFile();
            }
            catch (System.Exception ex)
            {
                Console.WriteLine(ex.Message);
                Console.WriteLine(ex.StackTrace);
            }
        }

        public static void ComputeSequenceFROverlap(string fastaFile, bool addReverse, string csvFileOut)
        {
            try
            {
                FileStream fs;
                try
                {
                    fs = new FileStream(fastaFile, FileMode.Open, FileAccess.Read, FileShare.Read);
                }
                catch (System.Exception)
                {
                    fs = new FileStream(fastaFile, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
                }

                Dictionary<char,Dictionary<string, long>> DicOfSeq = new Dictionary<char,Dictionary<string, long>>(30);
                for (int i = 0; i < 26; i++)
                    DicOfSeq.Add((char)('A' + i), new Dictionary<string, long>());

                using (StreamReader sr = new StreamReader(fs))
                {
                    string line;
                    while ((line = sr.ReadLine()) != null)
                    {
                        if (!line.StartsWith(">"))
                        {
                            if (DicOfSeq[line[0]].ContainsKey(line))
                                DicOfSeq[line[0]][line] += 1;
                            else
                                DicOfSeq[line[0]].Add(line, 1);

                            if (addReverse)
                            {
                                string rev = Reverse(line);

                                if (DicOfSeq[line[0]].ContainsKey(rev))
                                    DicOfSeq[line[0]][rev] += 1;
                                else
                                    DicOfSeq[line[0]].Add(rev, 1);
                            }                                
                        }
                    }
                }
                fs.Close();

                Dictionary<long, long> DicOfNb = new Dictionary<long, long>();
                for (long i = 1; i <= 40; i++)
                    DicOfNb.Add(i, 0);

                for(int i = 0; i < 26; i++)
                    foreach (long val in DicOfSeq[(char)('A' + i)].Values)
                        if (DicOfNb.ContainsKey(val))
                            DicOfNb[val] += 1;
                        else
                            DicOfNb.Add(val, 1);

                vsCSVWriter writer = new vsCSVWriter(csvFileOut);
                foreach (long key in DicOfNb.Keys)
                    writer.AddLine(key + "," + DicOfNb[key]);
                writer.writeToFile();
            }
            catch (System.Exception ex)
            {
                Console.WriteLine(ex.Message);
                Console.WriteLine(ex.StackTrace);
            }
        }

        public static void AppendUbiPredToMascotReport(string csvUbiFile, string csvMascotFile, string csvFileOut)
        {
            vsCSV csvUbi = new vsCSV(csvUbiFile);
            vsCSV csvMascot = new vsCSV(csvMascotFile);
            vsCSVWriter writer = new vsCSVWriter(csvFileOut);
            try
            {
                foreach(string lineMascot in csvMascot.LINES_LIST)
                {
                    string strToAppend = "";
                    try
                    {
                        string[] mSplits = lineMascot.Split(vsCSV._Generic_Separator);
                        if (mSplits.Length >= 17 && !lineMascot.StartsWith("Search"))
                        {
                            int indexStart = int.Parse(mSplits[16]);
                            int indexStop = int.Parse(mSplits[17]);
                            foreach (string lineUbi in csvUbi.LINES_LIST)
                            {
                                string[] splits = lineUbi.Split(vsCSV._Generic_Separator);
                                int indexUbi = int.Parse(splits[2]);
                                if (splits[4].Contains(mSplits[3]) && indexUbi >= indexStart && indexUbi <= indexStop)
                                    strToAppend += "," + splits[2] + "," + splits[3];
                            }
                        }
                    }
                    catch (System.Exception ex)
                    {
                        Console.WriteLine(ex.Message);
                        Console.WriteLine(ex.StackTrace);
                    }
                    writer.AddLine(lineMascot + strToAppend);
                }
                writer.writeToFile();
            }
            catch (System.Exception ex)
            {
                Console.WriteLine(ex.Message);
                Console.WriteLine(ex.StackTrace);
            }
        }

        public static void AppendhCKSAAPToMascotReport(string txtHCKSAAPFile, string csvMascotFile, string csvFileOut)
        {
            vsCSV csvUbi = new vsCSV(txtHCKSAAPFile);
            vsCSV csvMascot = new vsCSV(csvMascotFile);
            vsCSVWriter writer = new vsCSVWriter(csvFileOut);
            try
            {
                foreach (string lineMascot in csvMascot.LINES_LIST)
                {
                    string strToAppend = "";
                    try
                    {
                        string[] mSplits = lineMascot.Split(vsCSV._Generic_Separator);
                        if (mSplits.Length >= 17 && !lineMascot.StartsWith("Search"))
                        {
                            int indexStart = int.Parse(mSplits[16]);
                            int indexStop = int.Parse(mSplits[17]);
                            bool inIPI = false;
                            foreach (string lineUbi in csvUbi.LINES_LIST)
                            {
                                if(lineUbi.StartsWith(">"))
                                {
                                    if(inIPI)
                                        break;

                                    inIPI = false;
                                    if(lineUbi.StartsWith(">IPI:") && lineUbi.Contains(mSplits[3]))
                                        inIPI = true;
                                }
                                if (inIPI)
                                {
                                    string[] splits = lineUbi.Split('\t');
                                    int indexPos = -1;
                                    if (splits.Length > 6 && int.TryParse(splits[0], out indexPos))
                                    {
                                        if (indexPos >= indexStart && indexPos <= indexStop)
                                            strToAppend += "," + splits[0] + "," + splits[5];
                                    }
                                }
                            }
                        }
                    }
                    catch (System.Exception ex)
                    {
                        Console.WriteLine(ex.Message);
                        Console.WriteLine(ex.StackTrace);
                    }
                    writer.AddLine(lineMascot + strToAppend);
                }
                writer.writeToFile();
            }
            catch (System.Exception ex)
            {
                Console.WriteLine(ex.Message);
                Console.WriteLine(ex.StackTrace);
            }
        }
    }
}
