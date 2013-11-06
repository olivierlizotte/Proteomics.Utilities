using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteomics.Utilities
{
    public interface IConSol
    {
        void WriteLine(string line);
    }

    public class ConSolCommandLine : IConSol
    {
        public void WriteLine(string line)
        {
            Console.WriteLine(line);
        }
    }
}
