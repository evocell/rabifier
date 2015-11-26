# Rabifier

Rabifier is an automated bioinformatic pipeline for prediction and classification of Rab GTPases. 
For more detailed description of the pipeline check the references. 
If you prefer just to browse Rab GTPases in all sequenced Eukaryotic genomes visit [rabdb.org](http://rabdb.org).

Rabifier is freely distributed under the GNU General Public License, check the LICENCE file for details.

Please cite our papers if you use Rabifier in your projects.

* Rabifier2
* Thousands of Rab GTPases for the Cell Biologist. Diekmann Y, et al. PLoS Comput Biol 7(10): e1002217. 
[doi:10.1371/journal.pcbi.1002217](http://dx.plos.org/10.1371/journal.pcbi.1002217)

## Installation

Clone the git repository (including the rabifier-data submodule)

    git clone --recursive https://github.com/evocell/rabifier.git

### Python requirements, third party packages and other dependencies

Currently Rabifier supports Python 2.7.6, we are working on Python 3.4 support. 
Rabifier was tested only on GNU/Linux operating systems.

Rabifier requires additional third-party Python packages:

* biopython (1.65)
* scipy (0.16.0)
* matplotlib (1.4.3) (required for the seed database building)

Rabifier uses several bioinformatic tools, which are required for most of the classification stages. 
Create symlinks (or copy binaries) to the following programs to the 'bin' directory, or ensure that they are available in the system path.

* HMMER (3.1b1)
  * phmmer
  * hmmbuild
  * hmmpress
  * hmmscan
* BLAST+ (2.2.30)
  * blastp
* MEME4 (4.8.0)
  * meme
  * mast
* SUPFAM (1.75)
  * hmmscan.pl
  * ass3.pl
  * superfamily (NOTE: this is a folder containing several Superfamily database files)
* cd-hit (v4.6.1-2012-08-27)
  * cd-hit
* PRANK (v.140603)
  * prank
  
**NOTE: cd-hit and prank are only needed to build the seed database, they are not required for a regular Rabifier
usage with the precomputed database**

To install Superfamily models follow the instructions below, based on the [Superfamily website](http://supfam.org/SUPERFAMILY/howto_use_models.html).

    # Register at the Superfamily website to get username and password

    # Download files
    wget --http-user USERNAME --http-password PASSWORD -r -np -nd -e robots=off \
        -R 'index.html*' 'http://supfam.org/SUPERFAMILY/downloads/license/supfam-local-1.75/'
    wget http://scop.mrc-lmb.cam.ac.uk/scop/parse/dir.cla.scop.txt_1.75 -O dir.cla.scop.txt

    # Uncompress files
    gzip -d *.gz
    mv hmmlib_1.75 hmmblib

    # Make Perl scripts executable
    chmod u+x *.pl
    
    # Build the HMM library
    hmmpress hmmlib

### Seed database

Rabifier requires a seed database for Rab classification. A precomputed database is automatically cloned with a submodule. 
You can also create the database using rabifier/core.py on the raw manually curated datasets (also provided in the submodule).

## Usage

To predict Rabs in your sequences, save them in the FASTA format and then run:

    ./rabmyfire.py sequences.fa
    
For more options controlling Rabifier type:

    ./rabmyfire.py -h

## Bug reports and contributing

Please use the [issue tracker](https://github.com/evocell/rabifier/issues) to report bugs and suggest improvements.
