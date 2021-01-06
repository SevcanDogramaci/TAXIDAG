# TAXIDAG 

A comparative tool for taxon identification of low coverage ancient genomes.

TAXIDAG is an elaborative method  to distinguish between two closely related ancient genomes by utilizing mtDNA and finding different identical sites between these species.

TAXIDAG is developed and tested using Linux operating system (Ubuntu 18.04 LTS).

TAXIDAG is implemented with two similar species' ( Sheep and Goat ) genomes in this repository.

For further details, please look at the poster in docs folder.

### **Requirements**

→ **Working environment**

- Python 3+ (tested with 3.6)
- Bedtools ([https://bedtools.readthedocs.io/en/latest/](https://bedtools.readthedocs.io/en/latest/))

    Ubuntu installation : `sudo apt-get install bedtools`

→ **Required libraries**

### File Structure

**data/ :**

- Place reference files here !
- Place sample data in tmp folder within another folder with the name 1st_sample, 2nd_sample, 4th_sample etc.
- Place the final files produced by choar_first step. (trvposoar.bed, trvposchi.bed)

**docs/ :**

- Project documents available.

**choar_analysis.py :** 

- choar_analysis implementer file. Call this to run TAXIDAG.

**utils.py :** utility functions 

**settings.py :**

- File locations are set here !
- Change current_sample_dir to the directory name of sample data in tmp folder.
- Change aligned_sample_to_oar_file, aligned_sample_to_chi_file to sample data files' names.

Available in ```requirements.txt```. 
Install with ```pip install -r requirements.txt```.

### **How to install the program via GitHub**

```bash
git clone https://github.com/SevcanDogramaci/TAXIDAG.git
```

### How to run with debug mode

Use debug mode to see outputs of each step on the console. 

In addition to this, debug mode enables you to have each intermediate files be created on data_out folder. Otherwise, only needed files are generated.

In default debug mode is off. To activate debug mode run the following command :

```bash
python3 choar_anaylsis.py -d
```

