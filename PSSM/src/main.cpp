/*
 * Implement the PSSM for supervised motif finding.
 * Author: Yi Zhou
*/
#include "pssm.h"
using namespace std;

const float PSEUDOCOUNT = 0.25;

int main(int argc, char **argv)
{
    if (argc < 3 || argc > 4)
    {
        printf("[ERROR] %s takes 2 or 3 arguments, but %d were given.\n\n"
               "Usage: %s <motif_text_file> <DNA_fasta_file> [minimum_score_cutoff]\n\n"
               "This program implements the PSSM for supervised motif finding.\n"
               "If minimum_score_cutoff is omitted, max(0, lowest_score_among_sequences)\n"
               "in the training set is used.\n",
               argv[0], argc - 1, argv[0]);
        return 1;
    }
    // Read aligned motif sequences and the DNA sequence
    Svec motifs = preprocess::read_motifs(argv[1]);
    string DNA = preprocess::read_DNA(argv[2]);
    preprocess::normalize_sequence(DNA);
    for (auto &s : motifs)
    {
        preprocess::to_upper_case(s);
    }

    // Construct frequency matrix
    PSSM ans(motifs, DNA);
    printf("\nFrequency matrix:");
    ans.print_score_matrix(0);

    // Add pseudo-counts to the frequency matrix
    ans.add_pseudocount(PSEUDOCOUNT);

    // Convert to score matrix
    ans.convert_to_score_matrix();
    printf("\nPSSM:");
    ans.print_score_matrix(3);

    // Calculate scores for motif sequences in the input alignment
    printf("\nTraining set motif scores:\n");
    const size_t motif_num = motifs.size(),
                 motif_len = motifs[0].size();
    size_t min_motif_score_idx = 0;
    vector<float> input_motif_scores(motif_num);
    for (size_t i = 0; i < motif_num; i++)
    {
        input_motif_scores[i] = ans.calc_score_for_forward(motifs[i], 0);
        printf("\t%s\t%.3f\n", motifs[i].c_str(), input_motif_scores[i]);
        if (input_motif_scores[min_motif_score_idx] > input_motif_scores[i])
        {
            min_motif_score_idx = i;
        }
    }

    // Invert the PSSM to represent the complementary strand
    ans.generate_reverse_matrix();

    // Scan the direct and complementary strands for the motif
    const float MIN_SCORE = (argc == 4) ? stof(argv[3]) : input_motif_scores[min_motif_score_idx];

    printf("\nMatches with score %.3f or higher found in %s (length %zu bp):\n"
           "%-10s%-10s%-8s%-25s%-10s\n",
           MIN_SCORE, argv[2], DNA.size(),
           "Start", "End", "Strand", "Sequence", "Score");

    const size_t i_max = DNA.size() - motif_len + 1;
    for (size_t i = 0; i < i_max; i++)
    {
        float forward_score = ans.calc_score_for_forward(DNA, i),
              reverse_score = ans.calc_score_for_reverse(DNA, i);
        if (forward_score >= MIN_SCORE)
        {
            printf("%-10zu%-10zu%-8c%-25s%-10.3f\n",
                   i + 1, i + motif_len, '+',
                   DNA.substr(i, motif_len).c_str(),
                   forward_score);
        }
        if (reverse_score >= MIN_SCORE)
        {

            printf("%-10zu%-10zu%-8c%-25s%-10.3f\n",
                   i + 1, i + motif_len, '-',
                   ans.generate_reverse_strand(DNA.substr(i, motif_len)).c_str(),
                   reverse_score);
        }
    }
    return 0;
}
