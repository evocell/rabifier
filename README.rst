Rabifier
========

.. image:: https://img.shields.io/pypi/v/rabifier.svg
    :target: https://pypi.python.org/pypi/rabifier

Rabifier is an automated bioinformatic pipeline for prediction and classification of Rab GTPases. 
For more detailed description of the pipeline check the references. 
If you prefer just to browse Rab GTPases in all sequenced Eukaryotic genomes visit `rabdb.org <http://rabdb.org>`_.

Rabifier is freely distributed under the GNU General Public License, check the LICENCE file for details.

Please cite our papers if you use Rabifier in your projects.

* Rabifier2: an improved bioinformatic classifier of Rab GTPases. Surkont J, et al.
* Thousands of Rab GTPases for the Cell Biologist. Diekmann Y, et al. PLoS Comput Biol 7(10): e1002217.
  `doi:10.1371/journal.pcbi.1002217 <http://dx.plos.org/10.1371/journal.pcbi.1002217>`_

Installation
------------

To install Rabifier simply run

.. code-block:: bash

    pip install rabifier

Python requirements, third party packages and other dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rabifier supports Python 2.7 and Python 3.4.
Rabifier was tested only on a GNU/Linux operating system, we are not planning to support other platforms.

Rabifier depends on third-party Python libraries:

* biopython (>=1.66)
* numpy (>=1.10.1)
* scipy (>=0.16.1)

Rabifier uses several bioinformatic tools, which are required for most of the classification stages. 
Ensure that the following programs (or links pointing to them) are available in the system path.

* `HMMER <http://hmmer.janelia.org/>`_ (3.1b1): ``phmmer``, ``hmmbuild``, ``hmmpress``, ``hmmscan``
* `BLAST+ <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_ (2.2.30): ``blastp``
* `MEME4 <http://meme-suite.org/>`_ (4.10.2): ``meme``, ``mast``
* `Superfamily <http://supfam.cs.bris.ac.uk/SUPERFAMILY/>`_ (>=1.75): ``superfamily``
  (NOTE: this is a folder containing several Superfamily database files and scripts, see below)

If you have cloned this repository you need to compile the HMMs of Rab subfamilies using ``hmmpress``, i.e. run
``hmmpress rabifier/data/rab_subfamily.hmm``

Rabifier requires a seed database for Rab classification. A precomputed database is a part of this repository.
You can also create the database using ``rabifier-mkdb`` on the raw, manually curated data sets, available in
a seperate repository https://github.com/evocell/rabifier-data.
The build process requires additional software.

* `CD-HIT <http://weizhongli-lab.org/cd-hit/>`_ (v4.6.4): ``cd-hit``
* `PRANK <http://wasabiapp.org/software/prank/>`_ (v.150803): ``prank``
* `MAFFT <http://mafft.cbrc.jp/alignment/software/>`_ (v7.221): ``mafft``
* `matplotlib <http://matplotlib.org/>`_ (>=1.4.3) (optional)

To install Superfamily database follow the instructions below (based on the
`Superfamily website <http://supfam.org/SUPERFAMILY/howto_use_models.html>`_).

.. code-block:: bash

    # Register at the Superfamily website to get your username and password

    # Download files
    mkdir superfamily
    cd superfamily
    wget --http-user USERNAME --http-password PASSWORD -r -np -nd -e robots=off \
        -R 'index.html*' 'http://supfam.org/SUPERFAMILY/downloads/license/supfam-local-1.75/'
    wget http://scop.mrc-lmb.cam.ac.uk/scop/parse/dir.cla.scop.txt_1.75 -O dir.cla.scop.txt
    wget http://scop.mrc-lmb.cam.ac.uk/scop/parse/dir.des.scop.txt_1.75 -O dir.des.scop.txt

    # Uncompress files
    gzip -d *.gz
    mv hmmlib_1.75 hmmblib

    # Make Perl scripts executable
    chmod u+x *.pl
    
    # Build the HMM library
    hmmpress hmmlib

    # Create a symbolic link pointing to the database directory e.g. ln -s superfamily $HOME/bin/

Usage
-----

To run Rab prediction on protein sequences, save sequences in the
`FASTA format <https://en.wikipedia.org/wiki/FASTA_format>`_ and run:

.. code-block:: bash

    rabifier sequences.fa
    
For more options controlling Rabifier behaviour type:

.. code-block:: bash

    rabifier -h

Bug reports and contributing
----------------------------

Please use the `issue tracker <https://github.com/evocell/rabifier/issues>`_ to report bugs and suggest improvements.
