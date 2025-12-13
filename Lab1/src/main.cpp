#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <random>
#include <algorithm>
#include <unordered_set>
#include <time.h>
using namespace std;


class ColoringProblem
{
public:
    int GetRandom(int a, int b)
    {
        static mt19937 generator;
        uniform_int_distribution<int> uniform(a, b);
        return uniform(generator);
    }

    void ReadGraphFile(string filename)
    {
        ifstream fin(filename);
        string line;
        int vertices = 0, edges = 0;
        while (getline(fin, line))
        {
            if (line[0] == 'c')
            {
                continue;
            }

            stringstream line_input(line);
            char command;
            if (line[0] == 'p')
            {
                string type;
                line_input >> command >> type >> vertices >> edges;
                neighbour_sets.resize(vertices);
                colors.resize(vertices + 1);
            }
            else
            {
                int start, finish;
                line_input >> command >> start >> finish;
                // Edges in DIMACS file can be repeated, but it is not a problem for our sets
                neighbour_sets[start - 1].insert(finish - 1);
                neighbour_sets[finish - 1].insert(start - 1);
            }
        }
    }

    void GreedyGraphColoring()
    {
        //Number of vertices
        int n = neighbour_sets.size();

        //Colors array
        colors.assign(n, 0);

        //Current number of used colors
        maxcolor = 0;

        //Computing the degree of each vertex
        vector<int> deg(n);
        for (int i = 0; i < n; ++i)
            deg[i] = neighbour_sets[i].size();

        //Deletion order (reserving memory for n ints)
        vector<int> order;
        order.reserve(n);

        //Copying the graph for "deleting" vertices with bool markers
        vector<unordered_set<int>> adj_copy = neighbour_sets;
        vector<bool> removed(n, false);

        //n times for every vertex
        for (int step = 0; step < n; ++step)
        {
            //Finding the vertex with the smallest degree
            int v = -1;
            for (int i = 0; i < n; ++i)
            {
                if (!removed[i] && (v == -1 || adj_copy[i].size() < adj_copy[v].size()))
                    v = i;
            }

            //Removing vertex
            removed[v] = true;
            order.push_back(v);

            //Removing vertex from neighbors
            for (int u : adj_copy[v])
                adj_copy[u].erase(v);
            adj_copy[v].clear();
        }

        //Greedy coloring in reverse order
        reverse(order.begin(), order.end());
        for (int v : order)
        {
            //Collecting all colors used by neighbors
            unordered_set<int> used_colors;
            for (int u : neighbour_sets[v])
                if (colors[u] != 0)
                    used_colors.insert(colors[u]);

            //Finding the smallest available color and coloring the vertex
            int c = 1;
            while (used_colors.count(c)) ++c;
            colors[v] = c;
            if (c > maxcolor) maxcolor = c;
        }
    }


    bool Check()
    {
        for (size_t i = 0; i < neighbour_sets.size(); ++i)
        {
            if (colors[i] == 0)
            {
                cout << "Vertex " << i + 1 << " is not colored\n";
                return false;
            }
            for (int neighbour : neighbour_sets[i])
            {
                if (colors[neighbour] == colors[i])
                {
                    cout << "Neighbour vertices " << i + 1 << ", " << neighbour + 1 <<  " have the same color\n";
                    return false;
                }
            }
        }
        return true;
    }

    int GetNumberOfColors()
    {
        return maxcolor;
    }

    const vector<int>& GetColors()
    {
        return colors;
    }

private:
    vector<int> colors;
    int maxcolor = 1;
    vector<unordered_set<int>> neighbour_sets;
};

int main()
{
    vector<string> files = 
    { 
        "Graphs/myciel3.col", "Graphs/myciel7.col",
        "Graphs/school1.col", "Graphs/school1_nsh.col",
        "Graphs/anna.col", "Graphs/miles1000.col",
        "Graphs/miles1500.col",
        "Graphs/le450_5a.col", "Graphs/le450_15b.col",
        "Graphs/queen11_11.col"
    };
    ofstream fout("color.csv");
    fout << "Instance; Colors; Time (sec)\n";
    cout << "Instance; Colors; Time (sec)\n";
    for (string file : files)
    {
        ColoringProblem problem;
        problem.ReadGraphFile(file);
        clock_t start = clock();
        problem.GreedyGraphColoring();
        if (! problem.Check())
        {
            fout << "*** WARNING: incorrect coloring: ***\n";
            cout << "*** WARNING: incorrect coloring: ***\n";
        }
        fout << file << "; " << problem.GetNumberOfColors() << "; " << double(clock() - start) / 1000 << '\n';
        cout << file << "; " << problem.GetNumberOfColors() << "; " << double(clock() - start) / 1000 << '\n';
    }
    fout.close();
    return 0;
}