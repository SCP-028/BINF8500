#include <cstdio> // printf
#include "clust.h"
using namespace std;
const static size_t MAX_ITER = 2;

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
    const size_t NCOL = matrix[0].size() - 1;

    // Convert rows into objects for better readability
    vector<Sample> samples;
    for (size_t i = 1; i < NROW; i++) // skip header row
    {
        samples.emplace_back(matrix[i], i);
    }
    // Normalize data before clustering -> mean=0 and standard deviation=1 on each column
    kmeans::scale_features(samples);

    // Keep increasing k until the Bayesian Information Criterion is reached
    vector<float> BICs, WCSSs;
    // for (size_t k = 2; k < NROW; k++)
    for (size_t k = 2; k < 4; k++)
    {
        // Pick 1st centroid randomly
        srand(time(NULL)); // use time as seed
        int centroid_one = rand() % NROW;
        vector<Cluster> clusters;
        clusters.emplace_back(1);
        clusters[0].add_sample(samples[centroid_one]);
        // find the point that's the furthest from the last centroid and
        // assign it as the next centroid
        kmeans::initialize_clusters(clusters, samples, k);

        size_t iter_num = 0,
               num_of_changes = 1;
        while (iter_num < MAX_ITER && num_of_changes != 0)
        {
            num_of_changes = 0;
            // Assign each point to the closest centroid
            for (auto &sample : samples)
            {
                float min_distance = sample.get_distance_to_centroid();
                size_t current_cluster_id = sample.get_cluster_id(),
                       cluster_id_to_assign = current_cluster_id;
                vector<float> _features = sample.get_features();
                for (auto &cluster : clusters)
                {
                    vector<float> _centroid = cluster.get_centroid();
                    float _distance = kmeans::distance(_features, _centroid);
                    if (std::isnan(min_distance) || _distance < min_distance)
                    {
                        min_distance = _distance;
                        cluster_id_to_assign = cluster.get_cluster_id();
                    }
                }
                if (current_cluster_id != 0)
                {
                    if (cluster_id_to_assign != current_cluster_id)
                    {
                        clusters[current_cluster_id - 1].remove_sample(sample);
                        clusters[cluster_id_to_assign - 1].add_sample(sample);
                        num_of_changes++;
                    }
                }
                else
                {
                    clusters[cluster_id_to_assign - 1].add_sample(sample);
                    num_of_changes++;
                }
            }
            // Update the centroid by averaging all the points in the cluster
            for (auto &cluster : clusters)
            {
                if (cluster.get_cluster_size() == 0)
                {
                    size_t big_clust_idx = kmeans::cluster_with_max_size(clusters);
                    Sample furthest_sample = kmeans::furthest_sample_in_cluster(clusters[big_clust_idx]);
                    clusters[big_clust_idx].remove_sample(furthest_sample);
                    cluster.add_sample(furthest_sample);
                }
            }
            for (auto &cluster : clusters)
            {
                cluster.update_centroid();
            }
            iter_num++;
        }
        // Check the Bayesian Information Criterion
        float BIC, WCSS = 0.0;
        for (auto &cluster : clusters)
        {
            vector<Sample> ss = cluster.get_samples();
            for (auto &sample : ss)
            {
                WCSS += sample.get_distance_to_centroid();
            }
        }
        BIC = (log(NROW - 1) * k * NCOL) + WCSS; // BIC = ln(n) * kd + WCSS
        if (BICs.size() != 0 && BIC < BICs.back())
        {
            //////////////////////////////////////////////////////////////////////////////////////////
            printf("BIC turnpoint reached!");
            //////////////////////////////////////////////////////////////////////////////////////////
            break;
        }
        else
        {
            BICs.push_back(BIC);
            WCSSs.push_back(WCSS);
            for (auto &s : samples)
            {
                s.reset_sample();
            }
        }
    }
    for (size_t i = 0; i < BICs.size(); i++)
    {
        printf("k = %zu, WCSS = %.2f, BIC = %.2f\n\n", i + 2, WCSSs[i], BICs[i]);
    }
    return 0;
}
