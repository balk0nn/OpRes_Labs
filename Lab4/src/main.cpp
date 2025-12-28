#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <time.h>
#include <random>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include "../../Lab3/src/tabu.h"
using namespace std;


class BnBSolver {
public:
    void ReadGraphFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) throw std::runtime_error("Cannot open file");

        std::string line;
        int n = 0;

        // Чтение header p edge n m
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == 'c') continue;
            if (line[0] == 'p') {
                std::istringstream iss(line);
                std::string tmp;
                int m;
                iss >> tmp >> tmp >> n >> m;
                adjMatrix.assign(n, std::vector<bool>(n, false));
                break;
            }
        }

        // Чтение рёбер
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == 'c') continue;
            std::istringstream iss(line);
            char t;
            int u, v;
            iss >> t >> u >> v;
            if (t == 'e') {
                adjMatrix[u - 1][v - 1] = true;
                adjMatrix[v - 1][u - 1] = true;
            }
        }

        vertices.clear();
        for (int i = 0; i < n; ++i) vertices.push_back({i, 0});

        ClearAll();
    }

    void RunBnB()
    {
        // Current and best cliques
        Q.clear();
        Qmax.clear();


        // Color classes step count
        int n = vertices.size();
        C.assign(n + 1, std::vector<Vertex>());
        S.assign(n + 1, StepCount());
        for (auto &s : S) { s.i1 = 0; s.i2 = 0; }

        //depth
        level = 1;
        pk = 0;

        
        // Running initial heuristic to get a good lower bound
        MaxCliqueTabuSearch heuristic;
        vector<unordered_set<int>> neighbour_sets(n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (adjMatrix[i][j]) neighbour_sets[i].insert(j);
        heuristic.SetNeighbourSets(neighbour_sets);

        heuristic.RunSearch(100, 10);

        // Converting the result to BnBsolver format
        const auto &heurClique = heuristic.GetClique();
        Qmax.assign(heurClique.begin(), heurClique.end());

        // initial setup
        setDegrees(vertices);
        sortByDegree(vertices);
        initColors(vertices);

        // Starting the BnB
        BnBrecursion(vertices);
    }

    const std::vector<int>& GetClique() const { return Qmax; }

    bool Check() const {
        for (size_t i = 0; i < Qmax.size(); ++i)
            for (size_t j = i + 1; j < Qmax.size(); ++j)
                if (!adjMatrix[Qmax[i]][Qmax[j]]) return false;
        return true;
    }

    void ClearClique() {
        Q.clear();
        Qmax.clear();
        pk = 0;
        level = 1;
        for (auto &s : S) { s.i1 = 0; s.i2 = 0; }
        for (auto &cls : C) cls.clear();
    }

private:

    // Vertex structure for convenience
    struct Vertex {
        int id;
        int degree;
    };

    // for dynamic coloring if needed
    struct StepCount {
        int i1 = 0;
        int i2 = 0;
        void inc() { ++i1; }
    };

    std::vector<std::vector<bool>> adjMatrix;
    std::vector<Vertex> vertices;
    std::vector<int> Q, Qmax;
    std::vector<std::vector<Vertex>> C;
    std::vector<StepCount> S;
    int level = 0;
    int pk = 0;
    const float Tlimit = 0.025f;

    void ClearAll() {
        Q.clear();
        Qmax.clear();
        C.clear();
        S.clear();
        level = 1;
        pk = 0;
    }

    bool connection(int i, int j) const {
        return adjMatrix[i][j];
    }

    // counting the degrees of vertices
    void setDegrees(std::vector<Vertex>& R) {
        for (auto &v : R) {
            int d = 0;
            for (auto &u : R)
                if (connection(v.id, u.id)) ++d;
            v.degree = d;
        }
    }

    // sorting by degrees in descending order
    void sortByDegree(std::vector<Vertex>& R) {
        std::sort(R.begin(), R.end(), [](const Vertex &a, const Vertex &b) {
            return a.degree > b.degree;
        });
    }

    // Starting initial upper bound by coloring
    void initColors(std::vector<Vertex>& R) {
        int max_degree = 0;
        for (auto &v : R) if (v.degree > max_degree) max_degree = v.degree;
        for (size_t i = 0; i < R.size(); ++i) {
            if (i < (size_t)max_degree) R[i].degree = i + 1;
            else R[i].degree = max_degree + 1;
        }
    }

    // Checking if vertex v can be added to color class cls
    bool cut1(const Vertex &v, const std::vector<Vertex> &cls) const {
        for (auto &u : cls)
            if (connection(v.id, u.id)) return true;
        return false;
    }

    //Filtering vertices connected to the clique
    void cut2(const std::vector<Vertex> &A, std::vector<Vertex> &B) const {
        B.clear();
        const Vertex &last = A.back();
        for (size_t i = 0; i < A.size() - 1; ++i)
            if (connection(last.id, A[i].id)) B.push_back(A[i]);
    }

    // Greedy coloring + sorting for branch upper bound
    void color_sort(std::vector<Vertex> &R) {
        int j = 0;
        int maxno = 1;
        int min_k = Qmax.size() - Q.size() + 1;
        for (auto &cls : C) cls.clear();

        for (auto &v : R) {
            int k = 1;
            while (cut1(v, C[k])) k++;
            if (k > maxno) {
                maxno = k;
                if ((size_t)maxno + 1 >= C.size()) C.resize(maxno + 2);
            }
            C[k].push_back(v);
            if (k < min_k) R[j++] = v;
        }

        if (j > 0) R[j-1].degree = 0;
        if (min_k <= 0) min_k = 1;

        for (int k = min_k; k <= maxno; ++k)
            for (auto &v : C[k]) {
                R[j] = v;
                R[j++].degree = k;
            }
    }

    // The main BnB recursion function
    void BnBrecursion(std::vector<Vertex> R)
    {
        // Updating the depth statistic
        S[level].i1 += S[level-1].i1 - S[level].i2;
        S[level].i2 = S[level-1].i1;

        while (!R.empty()) {
            Vertex v = R.back();

            // UB check
            if ((int)Q.size() + v.degree > (int)Qmax.size()) {
                Q.push_back(v.id);

                std::vector<Vertex> Rp;
                cut2(R, Rp);

                if (!Rp.empty()) {

                    // Dynamic coloring condition
                    if ((float)S[level].i1 / ++pk < Tlimit) {
                        setDegrees(Rp);
                        sortByDegree(Rp);
                    }
                    color_sort(Rp);
                    S[level].inc();
                    level++;
                    BnBrecursion(Rp);
                    level--;

                // Found a leaf, check if better than best known
                } else if (Q.size() > Qmax.size()) {
                    Qmax = Q;
                }
                Q.pop_back();
                R.pop_back();

            // if not, nothing to search, go back
            } else {
                return;
            }
        }
    }
};

int main()
{
    //ios_base::sync_with_stdio(false);
    //cin.tie(nullptr);
    vector<string> files = 
    {
        "Graphs/brock200_1.clq",
        "Graphs/brock200_2.clq",
        "Graphs/brock200_3.clq",
        "Graphs/brock200_4.clq",
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
        "Graphs/p_hat1500-1.clq",
        "Graphs/p_hat300-3.clq",
        "Graphs/san1000.clq",
        "Graphs/sanr200_0.9.clq",
    };
    ofstream fout("clique_bnb.csv");
    fout << "File; Clique; Time (sec)\n";
    for (string file : files)
    {
        BnBSolver problem;
        problem.ReadGraphFile(file);
        problem.ClearClique();
        clock_t start = clock();
        problem.RunBnB();
        if (! problem.Check())
        {
            cout << "*** WARNING: incorrect clique ***\n";
            fout << "*** WARNING: incorrect clique ***\n";
        }
        fout << file << "; " << problem.GetClique().size() << "; " << double(clock() - start) / 1000 << '\n';
        cout << file << ", result - " << problem.GetClique().size() << ", time - " << double(clock() - start) / 1000 << '\n';
    }
    return 0;
}