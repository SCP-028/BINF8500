#include <cstdio> // std::printf
#include "clust.h"
using namespace std;
const static size_t MAX_ITER = 1000;

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

    // Convert rows into objects for better readability
    vector<Sample> samples;
    for (size_t i = 1; i < NROW; i++) // skip header row
    {
        Sample g(matrix[i], i);
        samples.push_back(g);
    }

    // Normalize data before clustering -> mean=0 and standard deviation=1 on each column
    kmeans::scale_features(samples);

    // Keep increasing k until the Bayesian Information Criterion is reached
    for (size_t k = 2; k < NROW; k++)
    {
        // Pick 1st centroid randomly
        srand(time(NULL)); // use time as seed
        int centroid_one = rand() % NROW;
        vector<Cluster> clusters;
        Cluster c(1, samples[centroid_one]);
        clusters.push_back(c);
        // find the point that's the furthest from the last centroid and
        // assign it as the next centroid
        for (size_t cluster_id = 2; cluster_id <= k; cluster_id++)
        {
            vector<float> last_centroid = clusters.back().get_centroid();
            float max_distance = 0.0;
            size_t new_centroid_sample_id;
            for (auto sample : samples)
            {
                if (sample.get_cluster_id() == 0) // doesn't have a cluster yet
                {
                    vector<float> _features = sample.get_features();
                    float _distance = kmeans::distance(_features, last_centroid);
                    if (_distance > max_distance)
                    {
                        max_distance = _distance;
                        new_centroid_sample_id = sample.get_sample_id();
                    }
                }
            }
            // create new cluster with this sample as the centroid
            Cluster c(cluster_id, samples[new_centroid_sample_id - 1]);
            clusters.push_back(c);
        }

        size_t iter_num = 0, num_of_changes;
        while (iter_num < MAX_ITER && num_of_changes != 0)
        {
            num_of_changes = 0;
            // Assign each point to the closest centroid
            for (auto sample : samples)
            {
                float min_distance = sample.get_distance_to_centroid();
                size_t current_cluster_id = sample.get_cluster_id(),
                       cluster_id_to_assign = 0;
                for (auto cluster : clusters)
                {
                    vector<float> _centroid = cluster.get_centroid(),
                                  _features = sample.get_features();
                    float _distance = kmeans::distance(_features, _centroid);
                    if (_distance < min_distance)
                    {
                        min_distance = _distance;
                        cluster_id_to_assign = cluster.get_cluster_id();
                    }
                }
                if (cluster_id_to_assign != current_cluster_id)
                {
                    clusters[current_cluster_id - 1].remove_sample(sample);
                    clusters[cluster_id_to_assign - 1].add_sample(sample);
                    num_of_changes++;
                }
            }
            // Update the centroid by averaging all the points in the cluster
            for (auto cluster : clusters)
            {
                if (cluster.get_samples().size() == 0)
                {
                    // deal with empty clusters
                }
                else
                {
                    cluster.update_centroid();
                }
            }
            // Check the Bayesian Information Criterion
        }
    }

    return 0;
}
