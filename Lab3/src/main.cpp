#include "tabu.h"


int main()
{
    int iterations;
    cout << "Number of iterations: ";
    cin >> iterations;
    int randomization;
    cout << "Randomization: ";
    cin >> randomization;
    
    vector<string> files = 
    { 
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
    
    ofstream fout("clique_tabu.csv");
    fout << "File; Clique; Time (sec)\n";
    
    for (string file : files)
    {
        MaxCliqueTabuSearch problem;
        problem.ReadGraphFile(file);
        clock_t start = clock();
        problem.RunSearch(iterations, randomization);
        
        if (!problem.Check())
        {
            cout << "*** WARNING: incorrect clique ***\n";
            fout << "*** WARNING: incorrect clique ***\n";
        }
        
        double time_taken = double(clock() - start) / CLOCKS_PER_SEC;
        fout << file << "; " << problem.GetClique().size() << "; " << time_taken << '\n';
        cout << file << ", result - " << problem.GetClique().size() << ", time - " << time_taken << " sec\n";
    }
    
    fout.close();
    return 0;
}
