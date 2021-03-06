---
Title: "Prenatal testosterone triggers long-term behavioral changes in male zebra finches: unravelling the neurogenomic mechanisms."
Authors: "Alex Bentz, Chad E. Niederhuth, Laura L. Carruth, and Kristen J. Navara"
Raw data: "GEO Accession: GSE100396"
---

This repository is for scripts and processed data for the paper:

Please cite this paper if you use any of the resources here.

All analyses performed on the Michigan State University High Performance Computing Cluster (HPCC)

To reproduce the analysis, follow these steps:

**NOTE #1:** This analysis assumes you will be using Anaconda and I have provided a yml file to easily create an environment for repeating analysies.

**1)** Clone this git repository

```
git clone https://github.com/niederhuth/Genomic-epigenomic-mechanisms-mediate-long-term-effects-of-prenatal-hormones-on-offspring-phenotype
cd Genomic-epigenomic-mechanisms-mediate-long-term-effects-of-prenatal-hormones-on-offspring-phenotype
```

**2)** Create the conda environment

```
conda env create -f scripts/zebra-finch-brains.yml
```

**3)** You will now need to create a symbolic link within this environment for methylpy to work. This will require you to cd into the environment located
in your anaconda (or miniconda) directory.

```
cd miniconda3/env/zebra-finch-brains/lib
ln -s libgsl.so.23.0.0 libgsl.so.0
```

**4)** Return to the cloned git repository

```
cd ~/Genomic-epigenomic-mechanisms-mediate-long-term-effects-of-prenatal-hormones-on-offspring-phenotype
```

**5)** Create a data folder and cd to it

```
mkdir data
cd data
```

**6)** Run the setup.sh script (see note #2)

```
bash ../scripts/setup.sh
```


