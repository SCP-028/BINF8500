#include <cstdio> // printf
#include "clust.h"
using namespace std;

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        printf("[Error] %s takes 1 arguement, but %d were given.\n\n"
               "Usage: %s <input_file>\n\n"
               "This program implements the k-means++ clustering algorithm.\n"
               "The appropriate value for k is chosen automatically using the\n"
               "Bayesian information criterion.\n",
               argv[0], argc - 1, argv[0]);
        return 1;
    }
    // Read file into vector of vectors
    vector<vector<string>> matrix;
    matrix.reserve(1000); // a little bit performance boost
    matrix = read_table(matrix, argv[1]);
    const size_t NROW = matrix.size();
    const size_t NCOL = matrix[0].size() - 1; // first column is used as rowname later

    // Convert rows into objects for better readability
    vector<Sample> samples;
    for (size_t i = 1; i < NROW; i++) // skip header row
    {
        Sample g(matrix[i]);
        samples.push_back(g);
    }

    // Normalize data before clustering -> mean 0 and sd 1 on each column
    scale_features(samples, NCOL);

    // Pick 1st centroid randomly, then pick the point that's furthest & repeat
    srand(time(NULL));                // Seed the time
    int centroid_one = rand() % NROW; // Generate the number, assign to variable.

    // Assign each point to the closest centroid

    // Update the location of the centroid by averaging all the points within the cluster

    // Check the Optimal Growth Temperature in the output (GC content & temp)
    for (auto g : samples)
    {
        printf("%s\t", g.sample_name.c_str());
        for (auto i : g.features)
        {
            printf("%.2f\t", i);
        }
        printf("\n");
    }
    return 0;
}
