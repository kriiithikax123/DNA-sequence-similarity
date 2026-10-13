# DNA Sequence Similarity in C++

A C++ implementation of generalized Jaccard similarity for DNA sequences using 2-bit nucleotide encoding, rolling k-mer encoding, and hash maps.

## Overview

This project calculates the similarity between two DNA sequences by comparing their k-mer frequency counts. A k-mer is a contiguous substring of length `k` in a DNA sequence.

The implementation uses 2-bit encoding to represent DNA nucleotides and rolling encoding to process consecutive k-mers efficiently. An `unordered_map` stores the frequency of each encoded k-mer.

The project also includes experiments to investigate how sequence length, mutations, and k-mer length affect similarity and execution time.

## Features

- Supports DNA bases `A`, `C`, `G`, and `T`.
- Encodes each nucleotide using 2 bits.
- Supports k-mer lengths from 1 to 32.
- Uses a 64-bit integer to represent encoded k-mers.
- Uses rolling encoding to update k-mers as the sequence is traversed.
- Counts k-mer occurrences using hash maps.
- Calculates generalized Jaccard similarity using k-mer frequencies.
- Generates random DNA sequences and introduces point mutations for experiments.
- Measures execution time and exports experiment results to CSV.

## Algorithm

### 1. K-mer counting

For a DNA sequence and a chosen k-mer length `k`, the program counts the occurrences of every contiguous k-mer.

For example, the sequence `ACGTAC` with `k = 3` contains:

- `ACG`
- `CGT`
- `GTA`
- `TAC`

Each k-mer is encoded as an integer, and its frequency is stored in a hash map.

### 2. Two-bit encoding

Each nucleotide is represented using two bits:

| Nucleotide | Encoding |
|---|---|
| A | `00` |
| C | `01` |
| G | `10` |
| T | `11` |

This allows a k-mer of up to 32 nucleotides to fit into a 64-bit unsigned integer.

### 3. Generalized Jaccard similarity

The similarity between two sequences is calculated from their k-mer frequency maps:

\[
J(A,B)=
\frac{\sum_x \min(c_A(x),c_B(x))}
{\sum_x \max(c_A(x),c_B(x))}
\]

Here, \(c_A(x)\) and \(c_B(x)\) represent the number of occurrences of k-mer \(x\) in sequences A and B.

K-mers that occur in only one sequence contribute to the denominator. Identical sequences have a similarity of 1 when the selected k-mer length is valid.

## Project Structure

```text
DNA-Sequence-Similarity/
├── src/
│   ├── main.cpp
│   ├── jaccard.cpp
│   └── jaccard.hpp
├── experiments/
├── results/
│   ├── similarity.csv
│   └── graphs/
├── README.md
└── .gitignore
```


## Requirements

- A C++ compiler that supports C++11 or later.
- A terminal or command prompt.
- GitHub is used to host the project.

The project was developed using C++ in Visual Studio Code on Windows with the MinGW-w64 GCC toolchain.

## Compilation and Execution

Open a terminal in the project root directory.

Compile the program:

```bash
g++ src/main.cpp src/jaccard.cpp -o main.exe
```

Run it on Windows PowerShell:

```powershell
.\main.exe
```

The program generates experiment results and writes them to a CSV file, according to the output path configured in the source code.

## Experiments

The experiments investigate three main factors:

1. **Sequence length:** How the execution time changes as DNA sequences become longer.
2. **Number of mutations:** How introducing point mutations affects the similarity between an original sequence and a mutated sequence.
3. **K-mer length:** How changing `k` affects similarity calculations and execution time.

The experiments use randomly generated DNA sequences, with mutations introduced incrementally. Execution time is measured for the Jaccard similarity calculation, while k-mer counting is performed outside the timed section.

The results are saved in CSV format and can be analyzed or visualized using spreadsheet software such as Microsoft Excel.

## Results and Visualization

The experiment results can be explored using the generated CSV file and graphs.

The visualizations cover:

- Sequence length versus runtime.
- Number of mutations versus Jaccard similarity.
- K-mer length versus runtime.

The measured results can be used to examine performance trends and the relationship between mutations and sequence similarity.

## Limitations

- The experiments use randomly generated DNA sequences rather than biological datasets.
- Runtime measurements can vary with hardware, compiler settings, and system load.
- Results depend on the chosen sequence lengths, k-mer lengths, and mutation levels.
- The similarity score compares k-mer frequency distributions and does not account for biological alignment or evolutionary relationships.

## Future Improvements

Potential extensions include:

- Adding automated unit tests.
- Comparing performance across different implementations.
- Testing on real DNA sequences.
- Improving experiment reproducibility through configurable random seeds.
- Adding automated graph generation.

## Author

Created as a C++ programming project exploring DNA sequence similarity, efficient k-mer encoding, and experimental performance analysis.
