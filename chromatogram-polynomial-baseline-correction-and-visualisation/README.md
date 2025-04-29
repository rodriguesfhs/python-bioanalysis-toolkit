# chromatogram_comparison

This script provides a function `plotter_chrom_compare` to load, baseline correct, and compare two chromatograms on the same plot. It is designed for analyzing data typically generated from chromatography experiments.

## Overview

The script defines a main function `plotter_chrom_compare` that orchestrates the following steps:

1.  **Data Loading:** Loads chromatogram data (`*_Chrom.csv`) and sample information (`Sample_ID.txt`) for two different experiments.
2.  **Data Slicing:** Allows you to specify a time range (`xmin`, `xmax`) to focus on a specific portion of each chromatogram.
3.  **Baseline Correction:** Applies a two-stage baseline correction: first a linear correction to address linear drift, followed by a polynomial correction to handle more complex baseline variations.
4.  **Plotting:** Generates a plot comparing the two baseline-corrected chromatograms, including labels, a legend, and a title derived from the sample information.

## File Structure

The script assumes a specific directory structure where the data files are located. For each experiment (A and B), it expects the following files within a subfolder named after the experiment number, inside a main folder identified by `folder`:


```
Data/CM/
└── [folderA]/
└── [exp_numberA]/
├── Sample_ID.txt
└── [exp_numberA]_Chrom.csv
└── [folderB]/
└── [exp_numberB]/
├── Sample_ID.txt
└── [exp_numberB]_Chrom.csv
```

**Note:** The current path `C:\Users\Documents\Data\CM\` is hardcoded within the `file_loader` function. **You will need to adjust this path** to match the location of your data on your system for the script to run correctly.

## `Sample_ID.txt` Format

The `Sample_ID.txt` file for each experiment is expected to be a colon-separated text file without a header. Each line contains a key-value pair (e.g., `Sample_Name:MySample`). The script extracts the sample name used for the plot label by assuming it is the second value on the second line of this file.

**Example `Sample_ID.txt`:**
```
Experiment:Control
Sample_Name:SampleX
Condition:Untreated
Analyst:Fabio
Date:2023-10-23
```

In this example, "SampleX" would be used as part of the plot label.

## `*_Chrom.csv` Format

The chromatogram data file (`[exp_number]_Chrom.csv`) is expected to be a space-delimited CSV file with two columns:

* **RT:** Retention Time (in some unit, later converted to minutes).
* **int:** Intensity.

## Usage

To use the `plotter_chrom_compare` function, you need to provide the following parameters for both chromatograms you want to compare:

```python
plotter_chrom_compare(
    folderA, exp_numberA, xminA, xmaxA, degreeA, xlimA, ylimA, colourA,
    folderB, exp_numberB, xminB, xmaxB, degreeB, xlimB, ylimB, colourB
)
```

**Parameters:**

* `folderA` (str): Folder identifier for the first experiment.
* `exp_numberA` (str): Experiment number identifier for the first experiment.
* `xminA` (float): Minimum time value for slicing the first chromatogram.
* `xmaxA` (float): Maximum time value for slicing the first chromatogram.
* `degreeA` (int): Degree of the polynomial for baseline fitting of the first chromatogram.
* `xlimA` (tuple): X-axis limits for the plot related to the first chromatogram (e.g., `(0, 30)`).
* `ylimA` (tuple): Y-axis limits for the plot related to the first chromatogram (e.g., `(-1000, 10000)`).
* `colourA` (str): Color for plotting the first chromatogram (e.g., `'blue'`, `'#FF0000'`).
* `folderB` (str): Folder identifier for the second experiment.
* `exp_numberB` (str): Experiment number identifier for the second experiment.
* `xminB` (float): Minimum time value for slicing the second chromatogram.
* `xmaxB` (float): Maximum time value for slicing the second chromatogram.
* `degreeB` (int): Degree of the polynomial for baseline fitting of the second chromatogram.
* `xlimB` (tuple): X-axis limits for the plot related to the second chromatogram.
* `ylimB` (tuple): Y-axis limits for the plot related to the second chromatogram.
* `colourB` (str): Color for plotting the second chromatogram.

Example:
```
plotter_chrom_compare(
    '181023', '430', 10, 25, 10, (10, 25), (-1000, 10000), 'dodgerblue',
    '181023', '431', 10, 20, 10, (10, 20), (0, 10000), 'orange'
)
```

**Important:** Ensure that the `folder` and `exp_number` parameters correctly point to the directories containing your `Sample_ID.txt` and `*_Chrom.csv` files. You will likely need to modify the hardcoded path in the `file_loader` function

Dependencies
This script relies on the following Python libraries:

* `pandas` for data manipulation.
* `matplotlib.pyplot` for plotting.
* `seaborn` for enhanced plot styling.
* `scipy.optimize` for curve fitting (baseline correction).
* `numpy` for numerical operations (polynomial fitting).
Make sure these libraries are installed in your Python environment before running the script. You can install them using pip:

```
pip install pandas matplotlib seaborn scipy numpy
```
