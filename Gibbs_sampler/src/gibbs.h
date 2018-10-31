#pragma once

#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

const static double PSEUDOCOUNT = 0.125;

template <typename T>
using Matrix = vector<vector<T>>;

namespace infiles {
inline Matrix<string> read_fasta(const string &input_file) {
    std::ifstream fin(input_file);
    std::stringstream buffer;
    if (fin.is_open()) {
        buffer << fin.rdbuf();
    } else {
        printf("Cannot open file: %s\n", input_file.c_str());
        exit(1);
    }
    string line;
    getline(buffer, line);  // first sequence
    vector<string> seq{line, ""};
    Matrix<string> ans;
    ans.reserve(30);
    while (getline(buffer, line)) {
        if (line[0] == '>') {
            ans.emplace_back(seq);
            seq[0] = line;
            seq[1] = "";
        } else {
            seq[1] += line;
        }
    }
    ans.emplace_back(seq);
    ans.shrink_to_fit();
    return ans;
}
inline Matrix<int> convert_nt_to_num(const Matrix<string> &fasta) {
    // fasta[0][0] is the sequence name, fasta[0][1] is the entire sequence in
    // one string
    const size_t seq_num = fasta.size();
    Matrix<int> ans(seq_num);
    for (size_t i = 0; i < seq_num; i++) {
        const size_t seq_len = fasta[i][1].size();
        ans[i].resize(seq_len);
        for (size_t j = 0; j < seq_len; j++) {
            switch (fasta[i][1][j]) {
                case 'A':
                    ans[i][j] = 0;
                    break;
                case 'C':
                    ans[i][j] = 1;
                    break;
                case 'G':
                    ans[i][j] = 2;
                    break;
                case 'T':
                    ans[i][j] = 3;
                    break;

                default:
                    break;
            }
        }
    }
    return ans;
}
}  // namespace infiles

namespace pssm {
inline Matrix<double> init_freq_matrix(const double pc,
                                       const double motif_len) {
    // Creates a `motif_len` by 4 matrix filled with pseudo-counts
    // as initial values
    Matrix<double> ans(motif_len, vector<double>(4, pc));
    return ans;
}
inline void count_bg_freq(const vector<int> &s, vector<double> &bg_freq) {
    for (auto &c : s) {
        bg_freq[c] += 1;
    }
    double total_score = 0.0;
    for (auto &score : bg_freq) {
        total_score += score;
    }
    for (auto &freq : bg_freq) {
        freq /= total_score;
    }
}
inline Matrix<double> calc_bg_freqs(const Matrix<int> &fasta) {
    const size_t seq_num = fasta.size();
    Matrix<double> bg_freqs(seq_num, vector<double>(4, PSEUDOCOUNT));
    for (size_t i = 0; i < seq_num; i++) {
        pssm::count_bg_freq(fasta[i], bg_freqs[i]);
    }
    return bg_freqs;
}
inline void count_char_occurrence(const vector<int> &s, const size_t s_idx,
                                  const size_t motif_len,
                                  Matrix<double> &freq_matrix) {
    //     pos      A  C  G  T
    //      0       pc pc pc pc
    //      1       pc pc pc pc
    //     ...      ...........
    // motif_len-1  pc pc pc pc
    for (size_t i = 0; i < motif_len; i++) {
        freq_matrix[i][s[i + s_idx]] += 1;
    }
}
inline void convert_to_pssm(const vector<double> &bg_freq,
                            Matrix<double> &freq_matrix) {
    // not taking log because can't have negative scores
    for (auto &row : freq_matrix) {
        double total_freq = 0.0;
        for (auto &ele : row) {
            total_freq += ele;
        }
        for (int i = 0; i < 4; i++) {
            row[i] = log2(row[i] / (total_freq * bg_freq[i]));
        }
    }
}
inline double calc_score(const Matrix<double> &pssm, const vector<int> &s,
                         const size_t start_pos, const size_t motif_len) {
    double total_score = 0.0;
    for (size_t i = 0; i < motif_len; i++) {
        total_score += pssm[i][s[i + start_pos]];
    }
    return total_score;
}
inline Matrix<double> build_pssm(const Matrix<int> &fasta,
                                 const Matrix<double> &bg_freqs,
                                 const size_t motif_len, const size_t seq_i,
                                 vector<size_t> &motif_positions) {
    Matrix<double> score_matrix =
        pssm::init_freq_matrix(PSEUDOCOUNT, motif_len);
    const size_t fasta_size = fasta.size();
    for (size_t i = 0; i < fasta_size; i++) {
        // use the background frequency of the left-out sequence , and use other
        // motifs to build the PSSM
        if (i != seq_i) {
            pssm::count_char_occurrence(fasta[i], motif_positions[i], motif_len,
                                        score_matrix);
        }
    }
    pssm::convert_to_pssm(bg_freqs[seq_i], score_matrix);
    return score_matrix;
}
}  // namespace pssm

namespace gibbs {
inline vector<size_t> init_motif_positions(const Matrix<int> &fasta,
                                           const size_t motif_len,
                                           mt19937 &generator) {
    vector<size_t> motif_positions;
    motif_positions.reserve(fasta.size());
    for (auto &seq : fasta) {
        uniform_int_distribution<size_t> dist(0, seq.size() - motif_len);
        motif_positions.emplace_back(dist(generator));
    }
    return motif_positions;
}
inline vector<double> calc_scores_in_seq(const Matrix<int> &fasta,
                                         const size_t seq_i,
                                         const size_t motif_len,
                                         const Matrix<double> &score_matrix) {
    const size_t motif_num = fasta[seq_i].size() - motif_len;
    vector<double> motif_scores(motif_num, 0.0);
    for (size_t _pos = 0; _pos < motif_num; _pos++) {
        motif_scores[_pos] =
            pssm::calc_score(score_matrix, fasta[seq_i], _pos, motif_len);
    }
    return motif_scores;
}
inline size_t update_position(const Matrix<int> &fasta,
                              const Matrix<double> &bg_freqs,
                              const size_t motif_len,
                              vector<size_t> &motif_positions,
                              vector<double> &fasta_scores,
                              mt19937 &generator) {
    size_t num_changes = 0;
    const size_t seq_num = fasta.size();
    for (size_t i = 0; i < seq_num; i++) {
        // Use all words except seq[i] to build the PSSM
        Matrix<double> score_matrix =
            pssm::build_pssm(fasta, bg_freqs, motif_len, i, motif_positions);

        // Use PSSM to assign score to each position in seq 1
        vector<double> motif_scores =
            gibbs::calc_scores_in_seq(fasta, i, motif_len, score_matrix);
        for (auto &ele : motif_scores) {
            ele = (ele < 0) ? 0 : ele;
        }
        // Select new location in seq 1 with probability proportional
        // to the previous scores using a weighted probability distribution
        std::discrete_distribution<size_t> dist(motif_scores.begin(),
                                                motif_scores.end());
        size_t new_pos = dist(generator);
        if (new_pos != motif_positions[i]) {
            num_changes++;
            motif_positions[i] = new_pos;
            fasta_scores[i] =
                pssm::calc_score(score_matrix, fasta[i], new_pos, motif_len);
        }
    }
    return num_changes;
}
inline void update_final_score(const vector<double> &fasta_scores,
                               const vector<size_t> &motif_positions,
                               double &max_score,
                               vector<size_t> &final_positions) {
    double total_score = 0.0;
    for (auto &ele : fasta_scores) {
        total_score += ele;
    }
    if (total_score > max_score) {
        max_score = total_score;
        final_positions = motif_positions;
    }
}
inline void shift_left_right(const Matrix<int> &fasta,
                             const Matrix<double> &bg_freqs,
                             const size_t motif_len,
                             vector<size_t> &motif_positions,
                             vector<double> &fasta_scores) {
    const size_t seq_num = fasta.size();
    for (size_t i = 0; i < seq_num; i++) {
        // Shift left 1nt if not already at position 0
        if (motif_positions[i] > 0) {
            motif_positions[i]--;
        }
        Matrix<double> score_matrix =
            pssm::build_pssm(fasta, bg_freqs, motif_len, i, motif_positions);
        double new_score = pssm::calc_score(score_matrix, fasta[i],
                                            motif_positions[i], motif_len);
        if (new_score > fasta_scores[i]) {
            fasta_scores[i] = new_score;
        } else {
            motif_positions[i]++;
        }
        // Check right boundary
        if (motif_positions[i] + motif_len + 2 <= fasta[i].size()) {
            motif_positions[i]++;
        }
        score_matrix =
            pssm::build_pssm(fasta, bg_freqs, motif_len, i, motif_positions);
        new_score = pssm::calc_score(score_matrix, fasta[i], motif_positions[i],
                                     motif_len);
        if (new_score > fasta_scores[i]) {
            fasta_scores[i] = new_score;
        } else {
            motif_positions[i]--;
        }
    }
}
inline void end_left_right(const Matrix<int> &fasta,
                           const Matrix<double> &bg_freqs, size_t &motif_len,
                           vector<size_t> &motif_positions,
                           vector<double> &fasta_scores) {
    size_t tmp_motif_len = motif_len, max_motif_len = motif_len;
    vector<size_t> tmp_motif_pos = motif_positions,
                   max_motif_pos = motif_positions;
    double tmp_score = 0.0, max_score = 0.0;
    const size_t seq_num = fasta.size();
    vector<double> tmp_scores(seq_num, 0.0), max_scores = fasta_scores;
    for (size_t i = 0; i < seq_num; i++) {
        max_score += fasta_scores[i];  // first max score is the initial score
    }
    // Right cut 1nt
    if (tmp_motif_len > 1) {
        tmp_motif_len--;
        for (size_t i = 0; i < seq_num; i++) {
            Matrix<double> score_matrix = pssm::build_pssm(
                fasta, bg_freqs, tmp_motif_len, i, tmp_motif_pos);
            tmp_scores[i] = pssm::calc_score(score_matrix, fasta[i],
                                             tmp_motif_pos[i], tmp_motif_len);
            tmp_score += tmp_scores[i];
        }
        if (tmp_score > max_score) {
            max_score = tmp_score;
            max_scores = tmp_scores;
            max_motif_pos = tmp_motif_pos;
            max_motif_len = tmp_motif_len;
        }
        tmp_motif_len++;
    }

    // Left cut 1nt
    tmp_score = 0.0;
    if (tmp_motif_len > 1) {
        tmp_motif_len--;
        for (size_t i = 0; i < seq_num; i++) {
            tmp_motif_pos[i]++;
            Matrix<double> score_matrix = pssm::build_pssm(
                fasta, bg_freqs, tmp_motif_len, i, tmp_motif_pos);
            tmp_scores[i] = pssm::calc_score(score_matrix, fasta[i],
                                             tmp_motif_pos[i], tmp_motif_len);
            tmp_score += tmp_scores[i];
        }
        if (max_score < tmp_score) {
            max_score = tmp_score;
            max_scores = tmp_scores;
            max_motif_pos = tmp_motif_pos;
            max_motif_len = tmp_motif_len;
        }
        tmp_motif_len++;
        tmp_motif_pos = motif_positions;
    }

    // Cut 1nt on both ends
    tmp_score = 0.0;
    if (tmp_motif_len > 2) {
        tmp_motif_len -= 2;
        for (size_t i = 0; i < seq_num; i++) {
            tmp_motif_pos[i]++;
            Matrix<double> score_matrix = pssm::build_pssm(
                fasta, bg_freqs, tmp_motif_len, i, tmp_motif_pos);
            tmp_scores[i] = pssm::calc_score(score_matrix, fasta[i],
                                             tmp_motif_pos[i], tmp_motif_len);
            tmp_score += tmp_scores[i];
        }
        if (max_score < tmp_score) {
            max_score = tmp_score;
            max_scores = tmp_scores;
            max_motif_pos = tmp_motif_pos;
            max_motif_len = tmp_motif_len;
        }
        tmp_motif_len += 2;
        tmp_motif_pos = motif_positions;
    }

    // Extend 1nt on the left
    bool at_left_end = false;
    for (size_t i = 0; i < seq_num; i++) {
        // Check left boundary
        if (tmp_motif_pos[i] == 0) {
            at_left_end = true;
        }
    }
    if (at_left_end != true) {
        tmp_motif_len++;
        for (size_t i = 0; i < seq_num; i++) {
            tmp_motif_pos[i]--;
            Matrix<double> score_matrix = pssm::build_pssm(
                fasta, bg_freqs, tmp_motif_len, i, tmp_motif_pos);
            tmp_scores[i] = pssm::calc_score(score_matrix, fasta[i],
                                             tmp_motif_pos[i], tmp_motif_len);
            tmp_score += tmp_scores[i];
        }
        if (max_score < tmp_score) {
            max_score = tmp_score;
            max_scores = tmp_scores;
            max_motif_pos = tmp_motif_pos;
            max_motif_len = tmp_motif_len;
        }
        tmp_motif_len--;
        tmp_motif_pos = motif_positions;
    } else {
        at_left_end = true;
    }

    // Extend 1nt on the right
    bool at_right_end = false;
    for (size_t i = 0; i < seq_num; i++) {
        // Check left boundary
        if (tmp_motif_pos[i] + tmp_motif_len + 2 > fasta[i].size()) {
            at_right_end = true;
        }
    }
    if (at_right_end != true) {
        tmp_motif_len++;
        for (size_t i = 0; i < seq_num; i++) {
            Matrix<double> score_matrix = pssm::build_pssm(
                fasta, bg_freqs, tmp_motif_len, i, tmp_motif_pos);
            tmp_scores[i] = pssm::calc_score(score_matrix, fasta[i],
                                             tmp_motif_pos[i], tmp_motif_len);
            tmp_score += tmp_scores[i];
        }
        if (max_score < tmp_score) {
            max_score = tmp_score;
            max_scores = tmp_scores;
            max_motif_pos = tmp_motif_pos;
            max_motif_len = tmp_motif_len;
        }
        tmp_motif_len--;
    } else {
        at_right_end = true;
    }

    // Extend 1nt on both ends
    if (at_left_end != true && at_right_end != true) {
        tmp_motif_len++;
        for (size_t i = 0; i < seq_num; i++) {
            tmp_motif_pos[i]--;
            Matrix<double> score_matrix = pssm::build_pssm(
                fasta, bg_freqs, tmp_motif_len, i, tmp_motif_pos);
            tmp_scores[i] = pssm::calc_score(score_matrix, fasta[i],
                                             tmp_motif_pos[i], tmp_motif_len);
            tmp_score += tmp_scores[i];
        }
        if (max_score < tmp_score) {
            max_score = tmp_score;
            max_scores = tmp_scores;
            max_motif_pos = tmp_motif_pos;
            max_motif_len = tmp_motif_len;
        }
        tmp_motif_len--;
        tmp_motif_pos = motif_positions;
    }
    motif_len = max_motif_len;
    motif_positions = max_motif_pos;
    fasta_scores = max_scores;
}

inline void final_scan(const Matrix<int> &fasta, const Matrix<double> &bg_freqs,
                       const size_t motif_len, vector<size_t> &motif_positions,
                       vector<double> &fasta_scores) {
    const size_t seq_num = fasta.size();
    for (size_t i = 0; i < seq_num; i++) {
        // Use all words except seq[i] to build the PSSM
        Matrix<double> score_matrix =
            pssm::build_pssm(fasta, bg_freqs, motif_len, i, motif_positions);

        // Use PSSM to assign score to each position in seq 1
        vector<double> motif_scores =
            gibbs::calc_scores_in_seq(fasta, i, motif_len, score_matrix);
        double high_score = 0.0;
        for (size_t pos = 0; pos < motif_scores.size(); pos++) {
            if (motif_scores[pos] > high_score) {
                high_score = motif_scores[pos];
                motif_positions[i] = pos;
            }
        }
        fasta_scores[i] = pssm::calc_score(score_matrix, fasta[i],
                                           motif_positions[i], motif_len);
    }
}
}  // namespace gibbs
