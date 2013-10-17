using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteomics.Utilities
{
    public static class MascotReportAnalyser
    {
        public static void MergeUbiSites(string csvFileA, string csvFileB, string csvOut)
        {
            int[] ubiPos = new int[]{34,36,38,40,42,44};
            vsCSV csvA = new vsCSV(csvFileA);
            vsCSV csvB = new vsCSV(csvFileB);
            vsCSVWriter writer = new vsCSVWriter(csvOut);
            int posUbiA = -1;
            int posUbiB = -1;
            foreach (string lineA in csvA.LINES_LIST)
            {
                string[] splitsA = lineA.Split(vsCSV._Generic_Separator);
                foreach(int indexUbiPos in ubiPos)
                {
                    if(splitsA.Length > indexUbiPos && int.TryParse(splitsA[indexUbiPos], out posUbiA))
                    {
                        foreach (string lineB in csvB.LINES_LIST)
                        {
                            string[] splitsB = lineB.Split(vsCSV._Generic_Separator);
                            
                            if(splitsB.Length > indexUbiPos && int.TryParse(splitsB[indexUbiPos], out posUbiB))
                            {
                                if(splitsA[3].CompareTo(splitsB[3]) == 0 && posUbiA == posUbiB)
                                {
                                    writer.AddLine( splitsA[3] + "," + vsCSV.GetFileName_NoExtension(csvFileA) + "," + splitsA[13] + "," + 
                                                        vsCSV.GetFileName_NoExtension(csvFileB) + "," + splitsB[13] + "," + 
                                                        posUbiA + "," + splitsA[indexUbiPos+1]);
                                }
                            }
                        }
                    }
                }
            }
            writer.writeToFile();
        }
    }
}
