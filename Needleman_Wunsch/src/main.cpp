/*
 * Implement the Needleman-Wunsch algorithm and align the two given sequences.
 * Author: Yi Zhou
*/
#include "nw.h"

using namespace std;

const unsigned PRINT_WIDTH = 80;

int main(int argc, char **argv)
{
    if (argc != 6)
    {
        printf("[Error] %s takes 5 arguments, but %d were given.\n\n"
               "Usage: %s <seq1> <seq2> <match_score> <mismatch_score> <gap_penalty>\n\n"
               "This program implements the Needleman-Wunsch algorithm and\n"
               "align the two given sequences.\n"
               "Matches and mismatches have uniform scores, whereas gaps have a linear penalty.\n",
               argv[0], argc - 1, argv[0]);
        return 1;
    }
    // Read file in and parse sequence into a single string
    const float MATCH = stof(argv[3]), MISMATCH = stof(argv[4]), GAP = stof(argv[5]);
    string seq1, seq2;
    seq1 = nw::read_fasta(argv[1]);
    seq2 = nw::read_fasta(argv[2]);
    nw::prepare_sequence(seq1);
    nw::prepare_sequence(seq2);

    // Construct score matrix (reverse the sequence)
    const size_t NROW = seq1.length(),
                 NCOL = seq2.length();
    vector<vector<float>> matrix(NROW + 1,
                                 vector<float>(NCOL + 1));
    nw::initialize_score_matrix(matrix, GAP);
    // Calculate score
    for (size_t i = 1; i <= NROW; i++)
    {
        for (size_t j = 1; j <= NCOL; j++)
        {
            matrix[i][j] = nw::max_score(nw::score_left(matrix, i, j, GAP),
                                         nw::score_top_left(matrix, i, j,
                                                            seq1, seq2,
                                                            MATCH, MISMATCH),
                                         nw::score_top(matrix, i, j, GAP));
        }
    }
    // Trace back and find the alignment

    // Print results
    nw::print_score_matrix(matrix, seq1, seq2);

    return 0;
}
