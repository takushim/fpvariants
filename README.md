# fpvariants

This repository is prepared as a supplemental material of our [review article for Frontiers in Physiology](https://www.frontiersin.org/journals/physiology).

## Usage
No additional procedures are necessary since all result files are output as pdf files in the `analyses` folder.

To rerun all procedures, install the following packages and run `run_all.ps1` after entering the `analyses` folder. It is highly recommended to install these packages in [a virtual environment of python](https://docs.python.org/3/library/venv.html).
* `numpy`
* `pandas`
* `matplotlib`

This repository contains two python scripts:
* `clcprep.py` - pre-processor for the csv files exported from [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
* `mapvariants.py` - script to plot variants on graphs.

## Author
[Takushi Miyoshi](https://github.com/takushim)

## License
Python scripts in this repository are licensed under the MIT licence. List of variants are exported from [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) as of Dec 3 2023.

