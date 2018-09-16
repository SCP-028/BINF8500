/*
 * Perform k-means++ clustering on a tab-delimited file with rows as samples
 * and columns as features. K values are automatically selected based on the
 * Bayesian Information Criterion.
 * Author: Yi Zhou
*/
#include <cstdio> // printf
#include <random>
#include "clust.h"
using namespace std;
const static size_t MAX_ITER = 1000, ITER_EACH = 100;

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
    random_device rn;   // obtain a random number
    mt19937 seed{rn()}; // seed the random number

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
    uniform_int_distribution<size_t> range(0, samples.size() - 1);

    // Normalize data before clustering -> mean=0 and standard deviation=1 on each column
    kmeans::scale_features(samples);

    // Keep increasing k until the turning point of the Bayesian Information Criterion is reached
    vector<float> BICs, WCSSs;
    vector<Cluster> ans;
    for (size_t k = 2; k < NROW; k++)
    {
        float _BIC = numeric_limits<float>::quiet_NaN(),
              _WCSS = numeric_limits<float>::quiet_NaN();
        vector<Cluster> _ans;
        /* repeat 10 times for each k and keep the one with min(BIC) because
         * different initial centroids may generate different clusters, and the
         * "turning point" of the BIC might not be the global optimum, so this
         * at least approaches the global optimum better
        */
        for (size_t i = 0; i < ITER_EACH; i++)
        {
            for (auto &s : samples)
            {
                s.reset_sample();
            }
            // Pick 1st centroid randomly
            vector<Cluster> clusters;
            clusters.emplace_back(1);
            clusters[0].add_sample(samples[range(seed)]);
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
                    float min_distance = std::numeric_limits<float>::max();
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
                        clusters[current_cluster_id - 1].remove_sample(sample);
                        clusters[cluster_id_to_assign - 1].add_sample(sample);
                        if (cluster_id_to_assign != current_cluster_id)
                        {
                            num_of_changes++;
                        }
                    }
                    else
                    {
                        clusters[cluster_id_to_assign - 1].add_sample(sample);
                        num_of_changes++;
                    }
                }
                // check for empty clusters
                Sample furthest_sample = kmeans::furthest_sample_in_clusters(clusters);
                for (auto &cluster : clusters)
                {
                    if (cluster.get_cluster_size() == 0)
                    {
                        clusters[furthest_sample.get_cluster_id() - 1].remove_sample(furthest_sample);
                        cluster.add_sample(furthest_sample);
                        furthest_sample = kmeans::furthest_sample_in_clusters(clusters);
                    }
                }
                // Update the centroid by averaging all the points in the cluster
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
            if (isnan(_BIC) || ((BIC < _BIC) && WCSS < _WCSS))
            {
                _BIC = BIC;
                _WCSS = WCSS;
                _ans = clusters;
            }
        }
        if (BICs.size() != 0 && _BIC > BICs.back())
        {
            printf("**************************\n"
                   "* BIC turnpoint reached! *\n"
                   "**************************\n"
                   "k = %zu, WCSS = %.2f, BIC = %.2f\n"
                   "Last BIC value was %.2f\n",
                   k, _WCSS, _BIC, BICs.back());
            break;
        }
        else
        {
            BICs.push_back(_BIC);
            WCSSs.push_back(_WCSS);
            ans = _ans;
            for (auto &s : samples)
            {
                s.reset_sample();
            }
        }
    }

    // print final result
    size_t final_k = BICs.size() + 1;
    printf("\n\nFinal result with k = %zu\n", final_k);
    for (auto &c : ans)
    {
        printf("\nCluster %zu\n----------\n", c.get_cluster_id());
        printf("\nSamples:\n\tDistance\tSample");
        for (auto &s : c.get_samples())
        {
            printf("\n\t%.2f\t%s",
                   s.get_distance_to_centroid(), s.get_sample_name().c_str());
        }
        printf("\n\n-----------------------------\n");
    }
    printf("\nCentroids\n----------\n");
    for (auto &c : ans)
    {
        printf("\nCluster %zu:", c.get_cluster_id());
        for (auto &p : c.get_centroid())
        {
            printf("\t%.2f", p);
        }
    }
    printf("\n\nMutual pairwise distances among centroids:\n\n");
    for (auto &c : ans)
    {
        printf("\tCluster %zu", c.get_cluster_id());
    }
    printf("\n");
    for (auto &c : ans)
    {
        for (size_t i = 0; i < ans.size(); i++)
        {
            vector<float> v1 = c.get_centroid(), v2 = ans[i].get_centroid();
            printf("\t%10.4f", kmeans::distance(v1, v2));
        }
        printf("\tCluster %zu\n", c.get_cluster_id());
    }
    printf("\n-----------------------------\n"
           "\nMean Distances Within Cluster\n"
           "-------------------------------\n");
    for (auto &c : ans)
    {
        float mean_distance = 0.0;
        for (auto &s : c.get_samples())
        {
            mean_distance += s.get_distance_to_centroid();
        }
        mean_distance /= c.get_samples().size();
        printf("\nCluster %zu\t%.2f", c.get_cluster_id(), mean_distance);
    }
    printf("\n\n-----------------------------\n\nk\tWCSS\tBIC");
    for (size_t i = 0; i < BICs.size(); i++)
    {
        printf("\n%zu\t%.2f\t%.2f", i + 2, WCSSs[i], BICs[i]);
    }
    return 0;
}
