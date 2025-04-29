# GC% Content Reducer (Randomized Codon Selection)

## Overview

This Python script, `GC_reducer_tool_rd.py`, is designed to optimize a given DNA sequence by reducing its Guanine-Cytosine (GC) content while ensuring that the translated protein sequence remains unchanged. It achieves this by employing a randomized codon selection strategy during the reverse translation process, favoring codons with lower GC content.

The script takes a DNA sequence as input, translates it into its corresponding protein sequence using a standard codon table, and then reverse translates this protein sequence back into a DNA sequence. During the reverse translation, for each amino acid, the script randomly selects one of its synonymous codons. This selection process inherently introduces variability in the resulting DNA sequence and, on average, leads to a reduction in GC content.

Finally, the script calculates and reports the initial and final GC percentages of the DNA sequence, the percentage reduction achieved, and performs a quality control check to confirm that the protein sequence is conserved after the optimization. The results are printed to the console and can be optionally saved as a CSV file and visualized as bar and line plots.

## Features

* **GC Content Reduction:** Optimizes DNA sequences to lower the overall GC percentage.
* **Protein Sequence Conservation:** Ensures that the amino acid sequence of the translated protein remains identical before and after optimization.
* **Randomized Codon Selection:** Uses a randomized approach when selecting synonymous codons for reverse translation.
* **Detailed Reporting:** Provides the initial and final GC percentages and the percentage of reduction.
* **Quality Control:** Verifies the conservation of the protein sequence.
* **Data Output:** Saves the results to a CSV file for further analysis (optional).
* **Visualizations:** Generates a bar chart comparing the initial and final GC percentages and a line plot showing the average GC content across the DNA sequence in 50-base pair windows (optional).

## Requirements

* Python 3.x
* `pandas` for creating and saving DataFrames.
* `numpy` (although imported, it might not be directly used in the final version of the provided code).
* `matplotlib` for generating plots.
* `seaborn` for enhanced plot aesthetics.
* `tqdm` for displaying progress bars during sequence processing.

You can install these dependencies using pip:

```bash
pip install pandas numpy matplotlib seaborn tqdm
```

## Usage
To use the GC_reducer_tool_rd.py script, you need to call the GC_reducer_tool_rd function with the following arguments:

```
GC_reducer_tool_rd(ID, DNA_seq, save_results_to, name_file, save_it)
```


### Arguments

* `ID` (str): A unique identifier for your DNA sequence (e.g., gene name). This will be used in the output file names and plot titles.
* `DNA_seq` (str): The DNA sequence you want to optimize. **Ensure it is provided in the 5' to 3' orientation**.
* `save_results_to` (str): The directory path where you want to save the output files (CSV and PNG images).
* `name_file` (str): The base name you want to use for the output files. The script will append suffixes like `_data.csv`, `_bar.png`, and `_GC_av.png` to this name.
* `save_it` (str): A string indicating whether to save the results. Use `'y'` or `'Y'` to save the output files. Any other value (e.g., `'n'` or `'N'`) will skip saving.

### Example
```
from GC_reducer_tool_rd import GC_reducer_tool_rd

ID = 'MyGene'
DNA_seq = 'ATGCGTAGCTAGCTAGCTAGCTAGCTAGCTAGC'
save_results_to = '/path/to/your/output/directory/'
name_file = 'my_gene_optimized'
save_it = 'y'

result_df = GC_reducer_tool_rd(ID, DNA_seq, save_results_to, name_file, save_it)

result_df
```

**Note:** Replace `/path/to/your/output/directory/` with the actual path where you want to save the results.

#### Output
When the script is executed, it will print the following information to the console:

* Sequence ID
* The optimized (GC%-reduced) DNA sequence
* The translated protein sequence
* The initial GC percentage
* The final GC percentage
* The percentage reduction in GC content
* Whether the protein sequence was conserved
If `save_it` is set to `'y'` or `'Y'`, the script will also generate the following files in the specified `save_results_to` directory:

* `{name_file}_data.csv`: A CSV file containing a transposed table of the results, including the sequence ID, optimized DNA, translated protein, initial and final GC percentages, GC reduction percentage, and protein sequence conservation status.
* `{name_file}_bar.png`: A bar chart comparing the initial and final GC percentages.
* `{name_file}_GC_av.png`: A line plot showing the average GC content in 50-base pair segments across the original and optimized DNA sequences.


Author: Dr Fabio H. S. Rodrigues  
Version: 13.0  
Date: 2025/04/29

