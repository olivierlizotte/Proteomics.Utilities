/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteomics.Utilities
{
    /// <summary>
    /// Max Flow optimization (Push relabel method)
    /// TODO Test the possibility of using this versus que QuickGraph library and set of Algorithms.
    /// </summary>
    public class PushRelabelMaximumFlow
    {
        int  NODES = 6;
        static int MIN(int X, int Y)
        {
            return X < Y ? X : Y;
        }

        static double INFINITE = 10000000;
 
        void push(List<List<int>> C, List<List<int>> F, int[] excess, int u, int v) {
                int send = MIN(excess[u], C[u][v] - F[u][v]);
                F[u][v] += send;
                F[v][u] -= send;
                excess[u] -= send;
                excess[v] += send;
        }
 
        void relabel(List<List<int>> C, List<List<int>> F, int[] height, int u) {
                int v;
                int min_height = int.MaxValue;
                for (v = 0; v < NODES; v++) {
                        if (C[u][v] - F[u][v] > 0) {
                                min_height = MIN(min_height, height[v]);
                                height[u] = min_height + 1;
                        }
                }
        }
 
        void discharge(List<List<int>> C, List<List<int>> F, int[] excess, int[] height, int[] seen, int u) {
                while (excess[u] > 0) {
                        if (seen[u] < NODES) {
                                int v = seen[u];
                                if ((C[u][v] - F[u][v] > 0) && (height[u] > height[v])){
                                        push(C, F, excess, u, v);
                                }
                                else
                                        seen[u] += 1;
                        } else {
                                relabel(C, F, height, u);
                                seen[u] = 0;
                        }
                }
        }
 
        void moveToFront(int i, int[] A) {
                int temp = A[i];
                int n;
                for (n = i; n > 0; n--){
                        A[n] = A[n-1];
                }
                A[0] = temp;
        }
 
        int pushRelabel(List<List<int>> C, List<List<int>> F, int source, int sink) 
        {
            int[] excess, height, list, seen;
            int i, p;
 
            excess = new int[6];
            height = new int[6];
            seen = new int[6];

            list = new int[6];
 
            for (i = 0, p = 0; i < NODES; i++){
                    if((i != source) && (i != sink)) {
                            list[p] = i;
                            p++;
                    }
            }
 
            height[source] = NODES;
            excess[source] = int.MaxValue;
            for (i = 0; i < NODES; i++)
                    push(C, F, excess, source, i);
 
            p = 0;
            while (p < NODES - 2) {
                    int u = list[p];
                    int old_height = height[u];
                    discharge(C, F, excess, height, seen, u);
                    if (height[u] > old_height) {
                            moveToFront(p,list);
                            p=0;
                    }
                    else
                            p += 1;
            }
            int maxflow = 0;
            for (i = 0; i < NODES; i++)
                    maxflow += F[source][i];
  
            return maxflow;
        }
 
        void printMatrix(List<List<int>> list) 
        {
            foreach(List<int> sublist in list)
            {
                string line = "";
                foreach(int item in sublist)
                    line += item + "\t";
                Console.WriteLine(line);
            }
        }
 
        public int test() 
        {
            List<List<int>> flow = new List<List<int>>(6);            
            List<List<int>> capacities = new List<List<int>>(6);

            for (int i = 0; i < 6; i++) 
            {
                List<int> newList = new List<int>();
                newList.Add(0);newList.Add(0);newList.Add(0);newList.Add(0);newList.Add(0);newList.Add(0);
                flow.Add(newList);
                newList = new List<int>();
                newList.Add(0); newList.Add(0); newList.Add(0); newList.Add(0); newList.Add(0); newList.Add(0);
                capacities.Add(newList);
            }
 
            //Sample graph
            capacities[0][1] = 2;
            capacities[0][2] = 9;
            capacities[1][2] = 1;
            capacities[1][3] = 0;
            capacities[1][4] = 0;
            capacities[2][4] = 7;
            capacities[3][5] = 7;
            capacities[4][5] = 4;
 
            Console.WriteLine("Capacity:");
            printMatrix(capacities);
 
            Console.WriteLine("Max Flow: " + pushRelabel(capacities, flow, 0, 5));
 
            Console.WriteLine("Flows:\n");
            printMatrix(flow);
 
            return 0;
        }//*/
    }
}
