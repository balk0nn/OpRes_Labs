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


class MaxCliqueTabuSearch
{
public:

    int GetRandom(int a, int b)
    {
        if (a >= b) return a;
        uniform_int_distribution<int> uniform(a, b);
        return uniform(rng);
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
                qco.resize(vertices);
                index.resize(vertices, -1);
                non_neighbours.resize(vertices);
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
        for (int i = 0; i < vertices; ++i)
        {
            for (int j = 0; j < vertices; ++j)
            {
                if (neighbour_sets[i].count(j) == 0 && i != j)
                    non_neighbours[i].insert(j);
            }
        }

        rng.seed((unsigned)time(nullptr));
    }

    void SetNeighbourSets(const std::vector<std::unordered_set<int>>& ns)
    {
        neighbour_sets = ns;

        int n = neighbour_sets.size();

        // Основные рабочие массивы
        qco.resize(n);
        index.assign(n, -1);
        tightness.assign(n, 0);

        // Строим множества несмежных вершин
        non_neighbours.resize(n);
        for (int i = 0; i < n; ++i) {
            non_neighbours[i].clear();
            for (int j = 0; j < n; ++j) {
                if (i != j && neighbour_sets[i].count(j) == 0)
                    non_neighbours[i].insert(j);
            }
        }

        // Сбрасываем текущее и лучшее решение
        best_clique.clear();
        q_border = 0;
        c_border = 0;

        // Инициализация qco: просто 0..n-1
        for (int i = 0; i < n; ++i) qco[i] = i;

        // Построение индексов
        RebuildIndex();

        // Инициализация tightness для пустой клики
        RebuildTightness();
    }

    void RunSearch(int starts, int randomization)
    {
        // Elite pool for path-relinking (keeps a few best solutions)
        vector<unordered_set<int>> elite;
        const int ELITE_MAX = 3;
        const int PR_FREQ = 5; // only attempt path-relinking every PR_FREQ starts

        // per-instance budgets
        // time_limit_seconds
        int n = (int)neighbour_sets.size();
        double time_limit_seconds = max(10.0, n * 0.05); // simple heuristic: 0.05s per vertex, min 10s
        if (time_limit_seconds > 120.0) time_limit_seconds = 120.0;
        clock_t global_start = clock();

        for (int iter = 0; iter < starts; ++iter)
        {
            // check global time budget
            double elapsed = double(clock() - global_start) / CLOCKS_PER_SEC;
            if (elapsed > time_limit_seconds) break;
            // Initialize working arrays
            ClearClique();
            for (size_t i = 0; i < neighbour_sets.size(); ++i)
            {
                qco[i] = i;
                index[i] = i;
            }

            // GRASP construction
            cur_randomization = randomization;
            RunInitialHeuristic(randomization);
            c_border = q_border;

            // Local search to improve the constructed solution (budget proportional to graph size)
            int ls_budget = max(20, n / 10);
            LocalSearch(ls_budget);

            // Capture current solution
            unordered_set<int> current = CaptureCurrentClique();

            // Path relinking: try combining with a random elite member occasionally
            if (!elite.empty() && (iter % PR_FREQ == 0))
            {
                int r = GetRandom(0, (int)elite.size() - 1);
                unordered_set<int> pr_candidate = PathRelink(current, elite[r]);
                if (pr_candidate.size() > current.size())
                {
                    current = pr_candidate;
                    RestoreClique(current);
                    RebuildIndex();
                    int ls_budget2 = max(10, n / 20);
                    LocalSearch(ls_budget2);
                }
            }

            // Shake / rebuild: perturb and re-search
            PerturbeClique(randomization);
            RebuildIndex();
            int ls_budget3 = max(10, n / 20);
            LocalSearch(ls_budget3);
            current = CaptureCurrentClique();

            // Insert into elite pool if good
            if (current.size() > 0)
            {
                // maintain unique elites
                bool inserted = false;
                for (auto &e : elite)
                {
                    if (e == current) { inserted = true; break; }
                }
                if (!inserted)
                {
                    elite.push_back(current);
                    sort(elite.begin(), elite.end(), [](const unordered_set<int>&a,const unordered_set<int>&b){return a.size()>b.size();});
                    if ((int)elite.size() > ELITE_MAX) elite.pop_back();
                }
            }

            // Update global best
            if (current.size() > best_clique.size())
            {
                best_clique = current;
            }
        }
    }

    unordered_set<int> CaptureCurrentClique()
    {
        unordered_set<int> s;
        for (int i = 0; i < q_border; ++i) s.insert(qco[i]);
        return s;
    }

    void RestoreClique(const unordered_set<int>& solution)
    {
        ClearClique();
        // place solution vertices at front of qco and fill rest with remaining vertices
        int n = (int)qco.size();
        int pos = 0;
        vector<char> in_sol(n, 0);
        for (int v : solution)
        {
            if (v < 0 || v >= n) continue;
            qco[pos] = v;
            index[v] = pos;
            in_sol[v] = 1;
            ++pos;
        }
        // fill remaining positions with vertices not in solution
        for (int v = 0; v < n; ++v)
        {
            if (in_sol[v]) continue;
            qco[pos] = v;
            index[v] = pos;
            ++pos;
        }
        // rebuild borders
        q_border = (int)solution.size();
        if (q_border > n) q_border = n;
        c_border = q_border;
        // ensure indices are consistent
        RebuildIndex();
        RebuildTightness();
    }

    // RebuildIndex: rebuild the `index` mapping from `qco` positions to
    // vertex ids. This is a cheap O(n) operation used after we rewrite
    // `qco` (e.g. in `RestoreClique`) to guarantee invariants used by
    // SwapVertices and other routines.
    void RebuildIndex()
    {
        int n = (int)qco.size();
        for (int i = 0; i < n; ++i)
        {
            int v = qco[i];
            if (v >= 0 && v < (int)index.size()) index[v] = i;
        }
    }

    void LocalSearch(int max_iterations)
    {
        int it = 0;
        int no_improve = 0;
        int best = q_border;
        while (it < max_iterations && no_improve < 25)
        {
            int before = q_border;
            // Greedy additions
            while (Move()) { if (q_border > best) { best = q_border; no_improve = 0; } }
            // Swaps
            if (!Swap1To1()) no_improve++; else { if (q_border>best) best=q_border; no_improve=0; }
            if (q_border == before) no_improve++;
            ++it;
        }
    }

    void PerturbeClique(int randomization)
    {
        if (q_border <= 1) return;
        int n = (int)neighbour_sets.size();
        double frac = (n > 1000) ? 0.07 : 0.12;
        int remove_cnt = max(1, (int)(q_border * frac));
        for (int r = 0; r < remove_cnt && q_border>1; ++r)
        {
            int idx = GetRandom(0, q_border - 1);
            int v = qco[idx];
            RemoveFromClique(v);
        }
        c_border = q_border;
        // try a brief rebuild
        RunInitialHeuristic(max(2, randomization/2));
        RebuildIndex();
        c_border = q_border;
    }

    unordered_set<int> PathRelink(const unordered_set<int>& A, const unordered_set<int>& B)
    {
        // Move from A towards B: iteratively try to add vertices from B\A when feasible
        unordered_set<int> current = A;
        unordered_set<int> target = B;
        unordered_set<int> best = current;
        bool improved = true;
        while (improved)
        {
            improved = false;
            // candidates in B not in current
            for (int v : target)
            {
                if (current.count(v)) continue;
                bool compatible = true;
                for (int u : current)
                    if (neighbour_sets[v].count(u) == 0) { compatible = false; break; }
                if (compatible)
                {
                    current.insert(v);
                    if (current.size() > best.size()) best = current;
                    improved = true;
                    break; // restart
                }
            }
        }
        return best;
    }

    const unordered_set<int>& GetClique()
    {
        return best_clique;
    }

    bool Check()
    {
        for (int i : best_clique)
        {
            for (int j : best_clique)
            {
                if (i != j && neighbour_sets[i].count(j) == 0)
                {
                    cout << "Returned subgraph is not clique\n";
                    return false;
                }
            }
        }
        return true;
    }

    void ClearClique()
    {
        q_border = 0;
        c_border = 0;
        // reset tightness cache
        tightness.assign(qco.size(), 0);
    }

private:
    vector<unordered_set<int>> neighbour_sets;
    vector<unordered_set<int>> non_neighbours;
    unordered_set<int> best_clique;
    vector<int> qco;
    vector<int> index;
    vector<int> tightness;
    mt19937 rng;

    // store last `randomization` parameter used;
    int cur_randomization = 1;
    int q_border = 0;
    int c_border = 0;

    int ComputeTightness(int vertex)
    {
        if (vertex >= 0 && vertex < (int)tightness.size()) return tightness[vertex];
        int t = 0;
        for (int i = 0; i < q_border; ++i)
        {
            if (non_neighbours[qco[i]].count(vertex))
                ++t;
        }
        return t;
    }

    // RebuildTightness: recompute the cached tightness (conflict counts)
    // for every vertex with respect to the current clique in `qco[0..q_border-1]`
    void RebuildTightness()
    {
        int n = (int)qco.size();
        tightness.assign(n, 0);
        for (int i = 0; i < q_border; ++i)
        {
            int u = qco[i];
            for (int v : non_neighbours[u])
                if (v >= 0 && v < n) tightness[v]++;
        }
    }

    void SwapVertices(int vertex, int border)
    {
        int n = (int)qco.size();
        if (border < 0 || border >= n) return;
        int pos_v = -1;
        if (vertex >= 0 && vertex < (int)index.size()) pos_v = index[vertex];
        if (pos_v < 0 || pos_v >= n)
        {
            // fallback: find vertex position in qco
            for (int k = 0; k < n; ++k)
            {
                if (qco[k] == vertex) { pos_v = k; index[vertex] = k; break; }
            }
            if (pos_v < 0 || pos_v >= n) return;
        }
        int vertex_at_border = qco[border];
        swap(qco[pos_v], qco[border]);
        // update indices
        index[qco[pos_v]] = pos_v;
        index[qco[border]] = border;
    }

    void InsertToClique(int i)
    {
        int n = (int)qco.size();
        int old_q = q_border;
        // for non-clique vertices, if they currently have zero conflicts, move them to candidate region
        for (int j : non_neighbours[i])
        {
            if (j < 0 || j >= n) continue;
            if (index[j] >= 0 && index[j] >= old_q)
            {
                if (tightness[j] == 0)
                {
                    --c_border;
                    SwapVertices(j, c_border);
                }
            }
        }
        SwapVertices(i, old_q);
        ++q_border;
        // rebuild tightness to be safe
        RebuildTightness();
    }

    void RemoveFromClique(int k)
    {
        int n = (int)qco.size();
        int old_q = q_border;
        for (int j : non_neighbours[k])
        {
            if (j < 0 || j >= n) continue;
            if (index[j] >= 0 && index[j] >= old_q)
            {
                if (tightness[j] == 1)
                {
                    SwapVertices(j, c_border);
                    c_border++;
                }
            }
        }
        --q_border;
        SwapVertices(k, q_border);
        // rebuild tightness to be safe
        RebuildTightness();
    }

    bool Swap1To1()
    {
        if (q_border <= 0) return false;
        // randomized scan order of clique vertices
        vector<int> order(q_border);
        for (int i = 0; i < q_border; ++i) order[i] = i;
        shuffle(order.begin(), order.end(), rng);
        for (int counter : order)
        {
            int vertex = qco[counter];
            // collect non-neighbors and randomize their order
            vector<int> candidates;
            for (int i : non_neighbours[vertex]) candidates.push_back(i);
            if (candidates.empty()) continue;
            shuffle(candidates.begin(), candidates.end(), rng);
            for (int i : candidates)
            {
                if (i >= 0 && i < (int)tightness.size() && tightness[i] == 1)
                {
                    RemoveFromClique(vertex);
                    InsertToClique(i);
                    return true;
                }
            }
        }
        return false;
    }

    bool Move()
    {
        if (c_border == q_border) return false;
        if (q_border < 0 || q_border >= (int)qco.size()) return false;
        // randomized selection among candidates in [q_border, c_border)
        int start = q_border;
        int end = min(c_border, (int)qco.size());
        if (start >= end) return false;
        vector<int> pos;
        pos.reserve(end - start);
        for (int i = start; i < end; ++i) pos.push_back(i);
        shuffle(pos.begin(), pos.end(), rng);
        for (int p : pos)
        {
            int vertex = qco[p];
            // pick vertices compatible (tightness==0)
            if (vertex >= 0 && vertex < (int)tightness.size() && tightness[vertex] == 0)
            {
                InsertToClique(vertex);
                return true;
            }
        }
        return false;
    }

    void RunInitialHeuristic(int randomization)
    {
        auto &generator = rng;
        vector<int> candidates(neighbour_sets.size());
        for (size_t i = 0; i < neighbour_sets.size(); ++i)
        {
            candidates[i] = i;
        }
        shuffle(candidates.begin(), candidates.end(), generator);
        while (! candidates.empty())
        {
            int last = candidates.size() - 1;
            int maxIdx;
            if (randomization <= 0) maxIdx = last;
            else maxIdx = min(randomization - 1, last);
            int rnd = GetRandom(0, maxIdx);
            int vertex = candidates[rnd];
            // Check compatibility with current clique (all vertices in qco[0..q_border-1])
            bool compatible = true;
            for (int i = 0; i < q_border; ++i)
            {
                if (neighbour_sets[vertex].count(qco[i]) == 0) { compatible = false; break; }
            }
            if (!compatible)
            {
                // remove this candidate and continue
                swap(candidates[rnd], candidates.back());
                candidates.pop_back();
                continue;
            }
            // accept vertex into the clique
            SwapVertices(vertex, q_border);
            ++q_border;
            // filter remaining candidates to those compatible with the whole clique
            for (int c = 0; c < (int)candidates.size(); ++c)
            {
                int candidate = candidates[c];
                bool ok = true;
                for (int i = 0; i < q_border; ++i)
                {
                    if (neighbour_sets[candidate].count(qco[i]) == 0) { ok = false; break; }
                }
                if (!ok)
                {
                    swap(candidates[c], candidates.back());
                    candidates.pop_back();
                    --c;
                }
            }
            shuffle(candidates.begin(), candidates.end(), generator);
        }
        // ensure tightness cache matches constructed clique
        RebuildTightness();
    }
};