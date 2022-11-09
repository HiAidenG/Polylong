# README
Sample data, utilities, and jupyter notebooks for the PolyLong project. 

More info below. 

# Setting up your environment

Create a new conda virtual environment:
```
conda create -n Polylong python=3.6
conda activate Polylong
```

Install DNABERT and dependencies:

```
conda install pytorch torchvision cudatoolkit=10.0 -c pytorch
git clone https://github.com/jerryji1993/DNABERT
cd DNABERT
python3 -m pip install --editable .
cd examples
python3 -m pip install -r requirements.txt
```

Install additional python packages:
```
conda install pandas scipy biopython seaborn matplotlib
```

You'll also need the seqkit and liftOver scripts. For those you have to add them to your PATH. \
liftOver is available here: http://hgdownload.cse.ucsc.edu/admin/exe/ \
seqkit ... ?

# WIP