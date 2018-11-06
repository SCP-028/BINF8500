/*
 * Implement the Gibbs sampler for unsupervised motif finding.
 * Author: Yi Zhou
 */
#include "gibbs.h"
using namespace std;

const static size_t INIT_SEED = 200, MAX_ITER = 2000, PLATEAU_CYCLES = 150;
// chance of: extending and shortening on the left side & on the right side
const static double PROB_SHIFT = 0.05, PROB_MOD_LEN = 0.10;

int main(int argc, char **argv) {
    if (argc != 3) {
        printf(
            "[ERROR] %s requires 2 arguments, but %d were given.\n"
            "Usage: %s <sequences.fasta> <estimated-motif-length>\n"
            "Example: %s data/E.coliRpoN-sequences-16-100nt.fasta 20\n",
            argv[0], argc - 1, argv[0], argv[0]);
        return 1;
    }
    random_device rn;
    // Read fasta file into matrix of strings, where:
    //     the first column is the name, and
    //     the second column is the sequence
    const Matrix<string> fasta_str = infiles::read_fasta(argv[1]);
    const Matrix<int> fasta = infiles::convert_nt_to_num(fasta_str);
    const Matrix<double> bg_freqs = pssm::calc_bg_freqs(fasta);

    // Global parameters
    uniform_real_distribution<double> shift_left_right(0, 1);

    const size_t fasta_size = fasta.size();
    vector<size_t> final_len(INIT_SEED);
    vector<float> max_scores(INIT_SEED);
    Matrix<size_t> final_positions(INIT_SEED);
#pragma omp parallel for
    for (size_t init_stat = 0; init_stat < INIT_SEED; init_stat++) {
        size_t motif_len = stoul(argv[2]);
        // If the guess is unreasonable, take its half as the starting point
        // The do...while loop takes care of motif lengths greater than the
        // sequence lengths
        while (motif_len >= fasta[0].size()) {
            motif_len /= 2;
        }

        double max_score = 0.0;
        vector<double> fasta_scores(fasta_size, 0.0);
        vector<size_t> final_position;

        // Initialization - randomly select a motif in each sequence
        mt19937 generator{rn()};
        vector<size_t> motif_positions =
            gibbs::init_motif_positions(fasta, motif_len, generator);

        size_t iter_num = 0, plateau_iter = 0;
        bool score_changed = false;
        do {
            iter_num++;
            // Update motif position in each sequence to a new position with
            // the scores as the probabilities
            gibbs::update_position(fasta, bg_freqs, motif_len, motif_positions,
                                   fasta_scores, generator);
            score_changed = gibbs::update_final_score(
                fasta_scores, motif_positions, max_score, final_position);

            // Shift each sequence left and right individually
            if (shift_left_right(generator) <= PROB_SHIFT) {
                gibbs::shift_left_right(fasta, bg_freqs, motif_len,
                                        motif_positions, fasta_scores);
                score_changed = gibbs::update_final_score(
                    fasta_scores, motif_positions, max_score, final_position);
            }
            // Shift all sequences
            if (shift_left_right(generator) <= PROB_MOD_LEN) {
                gibbs::end_left_right(fasta, bg_freqs, motif_len,
                                      motif_positions, fasta_scores);
                score_changed = gibbs::update_final_score(
                    fasta_scores, motif_positions, max_score, final_position);
            }
            if (score_changed == false) {
                plateau_iter++;
            } else {
                plateau_iter = 0;
            }
        } while (plateau_iter < PLATEAU_CYCLES && iter_num < MAX_ITER);
        // After MAX_ITER iterations, perform a final scan on the entire
        // sequence to guarantee the local maximum (a very important step)
        for (size_t i = 0; i < fasta_size; i++) {
            gibbs::final_scan(fasta, bg_freqs, motif_len, motif_positions,
                              fasta_scores);
        }
        gibbs::update_final_score(fasta_scores, motif_positions, max_score,
                                  final_position);
        // afterwards, check if shifting left / right would improve the score
        do {
            gibbs::end_left_right(fasta, bg_freqs, motif_len, motif_positions,
                                  fasta_scores);
            score_changed = gibbs::update_final_score(
                fasta_scores, motif_positions, max_score, final_position);
        } while (score_changed == true);
        // Store result of current seed
        max_scores[init_stat] = max_score;
        final_positions[init_stat] = final_position;
        final_len[init_stat] = motif_len;
    }

    // Find the max score in all different initial seeds
    double max_score = 0.0;
    vector<size_t> final_position;
    size_t motif_len = 0;
    for (size_t i = 0; i < INIT_SEED; i++) {
        if (max_scores[i] > max_score) {
            max_score = max_scores[i];
            final_position = final_positions[i];
            motif_len = final_len[i];
        }
    }
    if (max_score == 0.0) {
        printf("\nFailed to find a motif in all sequences.\n");
        return 1;
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
        const size_t motif_pos = final_position[i];
        const string motif_seq = fasta_str[i][1].substr(motif_pos, motif_len);
        printf("%s\t%zu-%zu\t%s\n", motif_seq.c_str(), motif_pos + 1,
               motif_pos + motif_len, fasta_str[i][0].c_str());
    }
    printf("\n\n");
    return 0;
}
