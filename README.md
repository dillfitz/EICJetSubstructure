# EIC jet substructure 

This repository is intended to setup a software workflow for using the EIC smear package to generate pythia truth MC and "smear" it using a fast detector simulation. 


### Requirements

This has only been tested on RCF, where access to the EIC environment is present. In order to use this, it is expected you have setup this environment via

```
setenv EIC_LEVEL pro 
source /afs/rhic.bnl.gov/eic/restructured/etc/eic_cshrc.csh

```

You need to also change the base path in the main file, `runSimWorkflow.py`, to point to the full path of wherever you checked this repo out. This is necessary so that python can find the relevant files.

### Usage
All that is necessary to use this is to run the main `runSimWorkflow.py` with the arguments as documented in the file. The output is a truth root file, which contains a tree with truth particle information from PYTHIA, and a smeared root file, which contains the particles after undergoing detector smearing.

The detector smearing is set in `smear/ePHENIXDetector.cpp`.

The output truth and smeared files can be befriended for further analysis within ROOT.


### Other notes

Don't be alarmed if you see "Error in <TROOT::TVector2>" type errors in the smearing. The smearing process sometimes removes a particle completely (e.g. for efficiency loss) and thus root complains when trying to smear a nonexistent particle.

### Example

To execute the runSimWorkflow.py, you can run the following:

python runSimWorkflow.py smeared truth 100 10 1 1000

This corresponds to running the workflow and generating a smeared file called “smeared.root”, a truth file called “truth.root”, proton energy = 100 GeV, electron energy = 10 Gev, minimum Q2 of 1 GeV^2,  and for 1000 events.