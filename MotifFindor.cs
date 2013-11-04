using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteomics.Utilities
{
    public static class MotifFindor
    {
        public static string[] LIR_Peptides = new string[]{ "SGGDDDWTHLSS",
                                                            "SASSEDYIIILP",
                                                            "LCGVSEWDPILE",
                                                            "ENEEDILVVTT",
                                                            "EGNSDMLVVTT",
                                                            "GSSEDSFVEIRM",
                                                            "AGLNSSWVELPM",
                                                            "ESLQGSWVELHF",
                                                            "ESDDDSYEVLDL",
                                                            "DSISGSWQAIQP",
                                                            "RVDHEEWEMVPR",
                                                            "ASSSFGWLSLDG",
                                                            "NEKALTWEEL",
                                                            "LSRPFTWEEI",
                                                            "SCDTDDFVMVPA",
                                                            "SCDTDDFVLVPH",
                                                            "HEDSDDFVLVPK",
                                                            "RSFEREYVVVEK",
                                                            "GNTHDDFVMIDF",
                                                            "DAHTFDFETIPH",
                                                            "DAATLTYDTLRF",
                                                            "LDGVGDWEDLQD",
                                                            "KIVDNDWLLPSY",
                                                            "VGYTPDWIFLLR",
                                                            "GSLEDDWDFLPP",
                                                            "EDEVDGWLIIDL",
                                                            "EKEDDEWILVDF",
                                                            "SPLLEDWDIISP",
                                                            "NSYRKEWEELFV",
                                                            "SSKDSGFTIVSP",
                                                            "PPDDAVFDIITD",
                                                            "EYRSRVYQMILE",
                                                            "EVRDRMWLKITI",
                                                            "LHPPSHWPLIKA"};
        /*
        public static int NbDEST(string sequence, int j)
        {
            if (sequence[j] == 'D' || sequence[j] == 'E' || sequence[j] == 'S' || sequence[j] == 'T')
                return 1;
            else
                return 0;
        }//*/

        public static int NbDEST(string sequence, int pos)
        {
            int nb = 0;
            for(int i = Math.Max(0, pos - 6); i < pos + 6 && i < sequence.Length; i++)
                if (sequence[i] == 'D' || sequence[i] == 'E' || sequence[i] == 'S' || sequence[i] == 'T')
                    nb++;
            return nb;
        }

        private static IEnumerable<string> GetAllLIRSequences(string DESTs, string WFY, string LIV, int length, string cumul)//, ref List<string> result)
        {
            if (cumul.Length == length)
                yield return cumul;//result.Add(cumul);
            else
            {
                if (cumul.Length == 6)
                {
                    GetAllLIRSequences(DESTs, WFY, LIV, length, cumul + "W");
                    GetAllLIRSequences(DESTs, WFY, LIV, length, cumul + "F");
                    GetAllLIRSequences(DESTs, WFY, LIV, length, cumul + "Y");
                }
                else if (cumul.Length == 9)
                {
                    GetAllLIRSequences(DESTs, WFY, LIV, length, cumul + "L");
                    GetAllLIRSequences(DESTs, WFY, LIV, length, cumul + "I");
                    GetAllLIRSequences(DESTs, WFY, LIV, length, cumul + "V");
                }
                else
                {
                    GetAllLIRSequences(DESTs, WFY, LIV, length, cumul + "D");
                    GetAllLIRSequences(DESTs, WFY, LIV, length, cumul + "E");
                    GetAllLIRSequences(DESTs, WFY, LIV, length, cumul + "S");
                    GetAllLIRSequences(DESTs, WFY, LIV, length, cumul + "T");
                    GetAllLIRSequences(DESTs, WFY, LIV, length, cumul + "x");
                }
            }
            //return result;
        }

        public static void ScoreSequences(string csvFileOut)
        {
            vsCSVWriter writer = new vsCSVWriter(csvFileOut);
            foreach (string motif in GetAllLIRSequences("DESTx", "WFY", "LIV", 12, ""))
            {
                double score = 0;
                foreach (string lir in LIR_Peptides)
                {
                    for (int i = 0; i < LIR_Peptides.Length; i++)
                        if (lir[i] == motif[i] || motif[i] == 'x')
                            score += 1;
                }
                if (score > 2)
                    writer.AddLine(motif + "," + score);                
            }
            writer.WriteToFile();
        }

        public static bool Score_LIR(string sequence, int i)
        {
            if (sequence[i - 2] == 'D' || sequence[i - 2] == 'E')
                if (sequence[i - 1] == 'D' || sequence[i - 1] == 'E' || sequence[i - 1] == 'S' || sequence[i - 1] == 'T')
                    if (sequence[i] == 'W' || sequence[i] == 'F' || sequence[i] == 'Y')
                        if (sequence[i + 1] == 'D' || sequence[i + 1] == 'E' || sequence[i + 1] == 'L' || sequence[i + 1] == 'I' || sequence[i + 1] == 'V')
                            if (sequence[i + 3] == 'L' || sequence[i + 3] == 'I' || sequence[i + 3] == 'V')
                                return true;
            return false;
        }

        public static double Score_LIRbkp(string sequence, int i)
        {
            double score = 0.0;
            int nbDest = NbDEST(sequence, i);
            if (sequence[i] == 'W' && nbDest > 0)
                score += 2;
            if (sequence[i] == 'F' && nbDest > 2)
            {
                score += 1;
                if (sequence[i + 1] == 'V' || sequence[i + 1] == 'C' || sequence[i + 1] == 'I' || sequence[i + 1] == 'E' || sequence[i + 1] == 'F')
                    score += 1;
            }
            if (sequence[i] == 'Y' && nbDest > 2)
                score += 1;
            if(score > 0)
            {
                if (sequence[i + 3] == 'L' || sequence[i + 3] == 'I' || sequence[i + 3] == 'V')
                {
                    score += 1;
                    score += nbDest * 0.2;
                    return score;
                }
            }
            return 0;
        }//*/

        public static void LIR(string fasta, string csvOut)
        {
            vsCSVWriter writer = new vsCSVWriter(csvOut);
            foreach (string[] protein in Fasta.FastaRead.GetSequences(fasta))
            //string header = "Known LIR";
            //foreach (string sequence in LIR_Peptides)
            {
                string header = protein[0];
                string sequence = protein[1];

                //if("W/F/Y" && pos+2 == "L/I/V" && "Enough E,D,S or T at +1 to -3")
                for (int i = 2; i + 3 < sequence.Length; i++)
                {
                    if(Score_LIR(sequence, i))
                    {
                    //double score = Score_LIR(sequence, i);
                    //int nbDest = NbDEST(sequence, i);
                    //if (score >= 2.6 && nbDest >= 3 && nbDest <= 7)
                    //{
                        writer.AddLine('"' + header + "\"," + i + "," + sequence.Substring(Math.Max(i - 9, 0), Math.Min(12, sequence.Length - Math.Max(i - 9, 0))) + "," + NbDEST(sequence, i));
                    }
                }
            }
            writer.WriteToFile();
        }
    }
}
