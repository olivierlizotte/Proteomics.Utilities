using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteomics.Utilities.Tools
{
    public static class AminoAcidTools
    {
        public static string Shuffle(string aa)
        {
            Random rd = new Random();
            List<char> listOfAA = new List<char>();
            foreach(char c in aa)
                listOfAA.Add(c);

            string shuffledAA = "";
            while (listOfAA.Count > 0)
            {
                int newIndex = rd.Next(0, listOfAA.Count - 1);
                shuffledAA += listOfAA[newIndex];
                listOfAA.RemoveAt(newIndex);
            } 
            return shuffledAA;
        }
    }
}
