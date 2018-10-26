/*
 * Implement the Gibbs sampler for unsupervised motif finding.
 * Author: Yi Zhou
 */
#include "gibbs.h"
using namespace std;

const static size_t MAX_ITER = 10000;

int main(int argc, char **argv) {
    if (argc != 3) {
        printf(
            "[ERROR] %s requires 2 arguments, but %d were given.\n"
            "Usage: %s <sequences.fasta> <motif-length>\n"
            "Example: %s data/E.coliRpoN-sequences-16-100nt.fasta 20\n",
            argv[0], argc - 1, argv[0], argv[0]);
        return 1;
    }
    // Random number generator
    random_device rn;
    mt19937 seed{rn()};

    // Read fasta file into matrix of strings, where:
    //     the first column is the name, and
    //     the second column is the sequence
    const Matrix<string> fasta = infiles::read_fasta(argv[1]);

    // Initialization - randomly select a motif in each sequence
    size_t motif_len = stoul(argv[2]);
    const size_t fasta_size = fasta.size();
    vector<size_t> motif_positions =
                       gibbs::init_motif_positions(fasta, motif_len, seed),
                   final_positions;

    vector<float> fasta_scores(fasta_size, 0.0);
    float max_score = 0.0;
    size_t num_changes = 1;
    size_t iter_num = 0;
    while (num_changes != 0 && iter_num <= MAX_ITER) {
        iter_num++;
        num_changes = gibbs::update_position(fasta, motif_len, seed,
                                             motif_positions, fasta_scores);
        gibbs::update_score(fasta_scores, motif_positions, max_score,
                            final_positions);
    }

    // Print results
    printf(
        "\nGibbs motif sampler output:\n\n"
        "\tInput file          : %s\n"
        "\tInitial motif length: %s\n"
        "\tFinal motif length  : %zu\n"
        "\tFinal score         : %.6f\n\n",
        argv[1], argv[2], motif_len, max_score);

    printf("Motif sequences and locations:\n\n");
    for (size_t i = 0; i < fasta_size; i++) {
        const size_t motif_pos = final_positions[i];
        const string motif_seq = fasta[i][1].substr(motif_pos, motif_len);
        printf("%s\t%zu-%zu\t%s\n", motif_seq.c_str(), motif_pos + 1,
               motif_pos + motif_len, fasta[i][0].c_str());
    }
    return 0;
}
