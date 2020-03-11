# StructuralErrorFinder

StructuralErrorFinder is a pipeline to find structural errors in your assembly based on your reads.

## Requirements
- Unix operating system (macOS or Linux)
- [Python 3](https://www.python.org)
- [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279671/) 

## Installation

The following instructions will install the latest version of StructuralErrorFinder:

```
git clone https://github.com/catrinehom/StructuralErrorFinder.git

cd StructuralErrorFinder/

chmod a+x StructuralErrorFinder.sh
chmod a+x PositionFinder.py
chmod a+x ErrorHandling.py
```

### Move to bin 
You might want to move the program to your bin to make the program globally excecutable. 
The placement of your bin depends on your system configuration, but common paths are:

```
/usr/local/bin/
```
OR
```
~/bin/
```

Example of move to bin:

```
mv StructuralErrorFinder.sh /usr/local/bin/
mv PositionFinder.py /usr/local/bin/
mv ErrorHandling.py /usr/local/bin/
```

## Usage

To run full pipeline:

```
./StructuralErrorFinder [-f <fastq reads>] [-a <fasta assembly] [-o <output name]
```

### Optional flags 
```
-m minimum alignment length for reads, default=5000
-l level of output (1, 2 or 3), default=1
```

