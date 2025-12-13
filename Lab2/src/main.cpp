#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <time.h>
#include <random>
#include <unordered_set>
#include <algorithm>
using namespace std;


class MaxCliqueProblem
{
public:
    static int GetRandom(int a, int b)
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

    void FindClique(int randomization, int iterations)
    {
        static mt19937 generator;
        for (int iteration = 0; iteration < iterations; ++iteration)
        {
            vector<int> clique;
            vector<int> candidates(neighbour_sets.size());
            for (int i = 0; i < neighbour_sets.size(); ++i)
            {
                candidates[i] = i;
            }
            shuffle(candidates.begin(), candidates.end(), generator);
            while (! candidates.empty())
            {
                int last = candidates.size() - 1;
                int rnd = GetRandom(0, min(randomization - 1, last));
                int vertex = candidates[rnd];
                clique.push_back(vertex);
                for (int c = 0; c < candidates.size(); ++c)
                {
                    int candidate = candidates[c];
                    if (neighbour_sets[vertex].count(candidate) == 0)
                    {
                        // Move the candidate to the end and pop it
                        swap(candidates[c], candidates[last]);
                        candidates.pop_back();
                        --c;
                    }
                }
                shuffle(candidates.begin(), candidates.end(), generator);
            }
            if (clique.size() > best_clique.size())
            {
                best_clique = clique;
            }
        }
    }

    const vector<int>& GetClique()
    {
        return best_clique;
    }

    bool Check()
    {
        if (unique(best_clique.begin(), best_clique.end()) != best_clique.end())
        {
            cout << "Duplicated vertices in the clique\n";
            return false;
        }
        for (int i : best_clique)
        {
            for (int j : best_clique)
            {
                if (i != j && neighbour_sets[i].count(j) == 0)
                {
                    cout << "Returned subgraph is not a clique\n";
                    return false;
                }
            }
        }
        return true;
    }

private:
    vector<unordered_set<int>> neighbour_sets;
    vector<int> best_clique;
};

int main()
{
    int iterations;
    cout << "Number of iterations: ";
    cin >> iterations;
    int randomization;
    cout << "Randomization: ";
    cin >> randomization;
    vector<string> files = {
        "Graphs/brock200_1.clq",
        "Graphs/brock200_2.clq",
        "Graphs/brock200_3.clq",
        "Graphs/brock200_4.clq",
        "Graphs/brock400_1.clq",
        "Graphs/brock400_2.clq",
        "Graphs/brock400_3.clq",
        "Graphs/brock400_4.clq",
        "Graphs/C125.9.clq",
        "Graphs/gen200_p0.9_44.clq",
        "Graphs/gen200_p0.9_55.clq",
        "Graphs/hamming8-4.clq",
        "Graphs/johnson16-2-4.clq",
        "Graphs/johnson8-2-4.clq",
        "Graphs/keller4.clq",
        "Graphs/MANN_a27.clq",
        "Graphs/MANN_a9.clq",
        "Graphs/p_hat1000-1.clq",
        "Graphs/p_hat1000-2.clq",
        "Graphs/p_hat1500-1.clq",
        "Graphs/p_hat300-3.clq",
        "Graphs/p_hat500-3.clq",
        "Graphs/san1000.clq",
        "Graphs/sanr200_0.9.clq",
        "Graphs/sanr400_0.7.clq"
    };
    ofstream fout("clique.csv");
    fout << "File; Clique; Time (sec)\n";
    for (string file : files)
    {
        MaxCliqueProblem problem;
        problem.ReadGraphFile(file);
        clock_t start = clock();
        problem.FindClique(randomization, iterations);
        if (! problem.Check())
        {
            cout << "*** WARNING: incorrect clique ***\n";
            fout << "*** WARNING: incorrect clique ***\n";
        }
        fout << file << "; " << problem.GetClique().size() << "; " << double(clock() - start) / 1000 << '\n';
        cout << file << ", result - " << problem.GetClique().size() << ", time - " << double(clock() - start) / 1000 << '\n';
    }
    fout.close();
    return 0;
}