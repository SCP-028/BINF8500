#include <cstdio>   // printf
#include <cstdlib>  // exit
#include <fstream>  // std::ifstream, ofstream
#include <iterator> // std::ostream_iterator
#include <string>   // std::string
#include <vector>   // std::vector, swap
using namespace std;
const size_t COLS = 4;

void swap_rows(vector<string> &v, size_t arow, size_t brow)
{
  /* Swap two rows in the vector. arow and brow are row numbers. */
  if (arow != brow)
  {
    for (unsigned i = 0; i < COLS; i++)
    {
      swap(v[arow * COLS + i], v[brow * COLS + i]);
    }
  }
}

void median_of_three(vector<string> &v, size_t start, size_t mid, size_t stop)
{
  /* Sort the first, middle and last elements in the vector, such that
  the smallest element is at the beginning, and the median is at the end. */
  if (v[start * COLS + 1] > v[mid * COLS + 1])
  {
    swap_rows(v, start, mid);
  }
  if (v[start * COLS + 1] > v[stop * COLS + 1])
  {
    swap_rows(v, start, stop);
  }
  if (v[stop * COLS + 1] > v[mid * COLS + 1])
  {
    swap_rows(v, mid, stop);
  }
}

size_t partition(vector<string> &v, size_t start, size_t stop)
{
  /* Perform the Hoare partition scheme. */
  size_t mid = (stop + start) / 2;
  median_of_three(v, start, mid, stop);
  string pvt = v[stop * COLS + 1];
  size_t i = start - 1, j = stop + 1;
  while (true)
  {
    do
    {
      i++;
    } while (v[i * COLS + 1] < pvt);

    do
    {
      j--;
    } while (v[j * COLS + 1] > pvt);

    if (i < j)
    {
      swap_rows(v, i, j);
    }
    else
    {
      return i;
    }
  }
}

void quicksort(vector<string> &v, size_t start, size_t stop)
{
  if (stop > start)
  {
    size_t pivot = partition(v, start, stop);

    if (pivot - start < stop - pivot)
    {
      quicksort(v, start, pivot - 1);
      quicksort(v, pivot, stop);
    }
    else
    {
      quicksort(v, pivot, stop);
      quicksort(v, start, pivot - 1);
    }
  }
}

int main(int argc, char **argv)
{
  if (argc == 1 || argc > 3)
  {
    printf("[ERR] Requires input FASTQ file and an optional output filename.\n\n"
           "Implements the QuickSort algorithm on the sequences inside a FASTQ file.\n\n"
           "Usage\n-----\n"
           "quicksort <input-fastq> [output-filename]\n");
    exit(1);
  }
  // Read file into vector
  vector<string> fastq;
  fastq.reserve(4000000);
  ifstream fin(argv[1]);
  if (fin.is_open())
  {
    string line;
    while (getline(fin, line))
    {
      fastq.push_back(line);
    }
    fastq.shrink_to_fit();
    fin.close();
  }
  else
  {
    printf("Cannot open file: %s\n", argv[1]);
    exit(1);
  }

  // Sort the vector `fastq`
  const size_t ROWS = fastq.size() / COLS;
  quicksort(fastq, 0, ROWS - 1);

  if (argc == 2)
  {

    for (size_t i = 0; i < fastq.size(); i++)
    {
      printf("%s\n", fastq[i].c_str());
    }
  }
  if (argc == 3)
  {
    // Write sorted vector into output file
    ofstream fout(argv[2]);
    if (fout.is_open())
    {
      ostream_iterator<string> iterout(fout, "\n");
      copy(fastq.begin(), fastq.end(), iterout);
      fout.close();
    }
    else
    {
      printf("Cannot write to file: %s\n", argv[2]);
      exit(1);
    }
  }
  return 0;
}
