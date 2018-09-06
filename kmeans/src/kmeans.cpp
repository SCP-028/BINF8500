#include <cstdio> // printf
#include "clust.h"
using namespace std;

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        printf("[Error] %s takes 1 argument, but %d were given.\n\n"
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
    matrix = kmeans::read_table(matrix, argv[1]);
    const size_t NROW = matrix.size();
    const size_t NCOL = matrix[0].size() - 1; // first column is used as row names later

    // Convert rows into objects for better readability
    vector<Sample> samples;
    for (size_t i = 1; i < NROW; i++) // skip header row
    {
        Sample g(matrix[i]);
        samples.push_back(g);
    }

    // Normalize data before clustering -> mean=0 and standard deviation=1 on each column
    kmeans::scale_features(samples);

    // Keep increasing k until the Bayesian information criterion is reached
    for (size_t k = 2; k < NROW; k++)
    {
        // Pick 1st centroid randomly
        srand(time(NULL)); // use time as seed
        int centroid_one = rand() % NROW;
        vector<Cluster> clusters;
        Cluster c(1, samples[centroid_one]);
        clusters.push_back(c);
        // pick the point that's furthest & repeat
        for (size_t cluster_id = 2; cluster_id <= k; cluster_id++)
        {
            // start with the only sample in the last cluster
            Sample furthest_sample = clusters[cluster_id - 1][0];
            for (auto item : samples)
            {
                if (item.get_cluster_id() == 0) // doesn't have a cluster yet
                {
                    float tmp_distance = kmeans::distance(item, clusters[0].get_centroid());
                    if (tmp_distance > furthest_sample.get_distance_to_centroid())
                    {
                        furthest_sample = item;
                    }
                }
            }
            Cluster c(cluster_id, furthest_sample);
            clusters.push_back(c);
        }
        // Assign each point to the closest centroid

        // Update the location of the centroid by averaging all the points within the cluster
    }

    // Check the Optimal Growth Temperature in the output (GC content & temp)
    for (auto g : samples)
    {
        printf("%s\t", g.get_sample_name().c_str());
        for (auto i : g.get_features())
        {
            printf("%.2f\t", i);
        }
        printf("\n");
    }
    return 0;
}
