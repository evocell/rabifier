#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part of Rabifier.
#
# Rabifier is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Rabifier is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rabifier.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function, division

import os
import argparse
import subprocess
import json
import math
import logging
import operator
import tempfile
import shutil
import csv

from Bio import SeqIO, pairwise2
import scipy.stats

from . import __version__
from . import config
from .utils import run_cmd, run_cmd_if_file_missing, Pathfinder


logger = logging.getLogger(__name__)


class Gprotein(object):

    def __init__(self, seqrecord, **kwargs):
        self.seqrecord = seqrecord
        self.bh_evalue = kwargs.get('bh_evalue', config['param']['bh_evalue'])
        self.motif_number = kwargs.get('motif_number', config['param']['motif_number'])

        self.gdomain_regions = []
        self.evalue_bh_rabs = None
        self.evalue_bh_non_rabs = None
        self.rabf_motifs = []

        self.rab_subfamily = None
        self.rab_subfamily_score = None
        self.rab_subfamily_top5 = []

    def __str__(self):
        return self.summarize()

    def is_gprotein(self):
        return True if self.gdomain_regions else False

    def is_best_hit_rab(self):
        if self.evalue_bh_rabs is not None:
            if self.evalue_bh_rabs < self.bh_evalue and \
                    (self.evalue_bh_non_rabs is None or self.evalue_bh_rabs <= self.evalue_bh_non_rabs):
                return True
        else:
            return False

    def has_rabf_motif(self):
        """Checks if the sequence has enough RabF motifs within the G domain

        If there exists more than one G domain in the sequence enough RabF motifs is required in at least one
        of those domains to classify the sequence as a Rab.
        """

        if self.rabf_motifs:
            for gdomain in self.gdomain_regions:
                beg, end = map(int, gdomain.split('-'))
                motifs = [x for x in self.rabf_motifs if x[1] >= beg and x[2] <= end]
                if motifs:
                    matches = int(pairwise2.align.globalxx('12345', ''.join(str(x[0]) for x in motifs))[0][2])
                    if matches >= self.motif_number:
                        return True
        return False

    def is_rab(self):
        return True if all([self.is_gprotein(), self.is_best_hit_rab(), self.has_rabf_motif()]) else False

    def to_dict(self):
        d = {
            'id': self.seqrecord.id,
            'sequence': {'header': self.seqrecord.description, 'seq': str(self.seqrecord.seq)},
            'gdomain_regions': self.gdomain_regions,
            'evalue_bh_rabs': self.evalue_bh_rabs,
            'evalue_bh_non_rabs': self.evalue_bh_non_rabs,
            'rabf_motifs': self.rabf_motifs,
            'is_rab': self.is_rab(),
            'rab_subfamily': [self.rab_subfamily, self.rab_subfamily_score],
            'rab_subfamily_top_5': self.rab_subfamily_top5
        }
        return d

    def to_list(self):
        return [
            self.seqrecord.id,
            ' '.join(self.gdomain_regions),
            self.evalue_bh_rabs,
            self.evalue_bh_non_rabs,
            'None' if not self.rabf_motifs else '|'.join(' '.join(map(str, x)) for x in self.rabf_motifs),
            self.is_rab(),
            self.rab_subfamily,
            self.rab_subfamily_score,
            '|'.join(' '.join(map(str, x)) for x in self.rab_subfamily_top5)
        ]

    def summarize(self):
        """ G protein annotation summary in a text format

        :return: A string summary of the annotation
        :rtype: str
        """
        data = [
            ['Sequence ID', self.seqrecord.id],
            ['G domain', ' '.join(self.gdomain_regions) if self.gdomain_regions else None],
            ['E-value vs rab db', self.evalue_bh_rabs],
            ['E-value vs non-rab db', self.evalue_bh_non_rabs],
            ['RabF motifs', ' '.join(map(str, self.rabf_motifs)) if self.rabf_motifs else None],
            ['Is Rab?', self.is_rab()]
        ]
        summary = ''
        for name, value in data:
            summary += '{:25s}{}\n'.format(name, value)
        if self.is_rab():
            summary += '{:25s}{}\n'.format('Top 5 subfamilies',
                                           ', '.join('{:s} ({:.2g})'.format(name, score) for name, score
                                                     in self.rab_subfamily_top5))
        return summary


class Phase1(object):
    """Determines which proteins, present in the fasta file, belong to the Rab family

    1) Determine if the G domain is present in the sequence.
    2) Check if the sequence's best hit against the reference Rab database is below the chosen e-value threshold and
    if the e-value is lower than the e-value of the best hit against the reference non-Rab database.
    3) Check if enough Rab motifs are present in the sequence

    :param str fastafile: Path to a fasta file with protein sequences
    :raises IOError: if fastafile is missing/empty or not in a valid fasta format
    :raises ValueError: if sequence is not a protein sequence
    :return: SeqRecord list
    """

    def __init__(self, fastafile, tmp, **kwargs):
        self.fastafile = fastafile

        # Seed database
        def get_db_file(name):
            path = os.path.join(kwargs.get('seed_path', config['seed']['path']), config['seed'][name])
            if os.path.exists(path):
                return path
            else:
                raise IOError("File '{}' doesn't exist".format(path))
        self.rabs_db = kwargs.get('rab_db', get_db_file('rab_db'))
        self.non_rabs_db = kwargs.get('non_rab_db', get_db_file('non_rab_db'))
        self.rabf_motifs = kwargs.get('rab_f', get_db_file('rab_f'))
        # TODO make a better database check for integrity and correctness
        for f in (self.rabs_db, self.non_rabs_db, self.rabf_motifs):
            if not os.path.exists(f):
                logger.error('Database file {} is missing.'.format(f))
                raise RuntimeError('Database file {} is missing.'.format(f))

        # Technical parameters
        self.fast = kwargs.get('fast', config['param']['fast'])
        self.cpu = str(kwargs.get('cpu', config['param']['cpu']))

        # Classification parameters
        self.bh_evalue = kwargs.get('bh_evalue', config['param']['bh_evalue'])
        self.motif_evalue = kwargs.get('motif_evalue', config['param']['motif_evalue'])

        # Set paths
        self.tmp = tmp
        self.tmpfname = os.path.join(tmp, os.path.basename(fastafile))  # base name for output files
        self.pathfinder = Pathfinder(True)
        self.pathfinder.add_path(self.pathfinder['superfamily'])

        self.query_id_list = [seq.id for seq in SeqIO.parse(self.fastafile, 'fasta')]
        self.gproteins = {}
        for x in SeqIO.parse(self.fastafile, 'fasta'):
            #  FIXME check if the sequence is a protein sequence
            self.gproteins[x.id] = Gprotein(
                x,
                bh_evalue=self.bh_evalue,
                motif_number=kwargs.get('motif_number', config['param']['motif_number'])
            )
        if not self.gproteins:
            raise IOError('IOError. File {} is missing/empty or has an incorrect format.\n'.format(fastafile))

    def __call__(self):
        # determine if sequences have G domains
        protein2gdomain = self.find_gdomain()
        for prot_id, gdomain in protein2gdomain.items():
            self.gproteins[prot_id].gdomain_regions = gdomain

        if self.fast:
            n = self._update_fasta(protein2gdomain.keys())
            if n == 0:
                return 0, ''

        # check how well sequences score against the rab and non-rab databases
        for prot_id, score in self.find_best_hit().items():
            self.gproteins[prot_id].evalue_bh_rabs = score['evalue_bh_rabs']
            self.gproteins[prot_id].evalue_bh_non_rabs = score['evalue_bh_non_rabs']

        # confirm presence of RabF-motifs
        for prot_id, motifs in self.find_rabf_motifs().items():
            self.gproteins[prot_id].rabf_motifs = motifs

        return 0, ''

    def find_gdomain(self):
        run_cmd([self.pathfinder['hmmscan.pl'], '-o', self.tmpfname + '.supfam_scan', '-E', '10', '-Z', '15438',
                 self.pathfinder['hmmlib'], self.fastafile, '--hmmscan', self.pathfinder['hmmscan'], '--threads',
                 self.cpu, '--tempdir', self.tmp])
        run_cmd([self.pathfinder.get('ass3.pl'), '-t', 'n', '-f', self.cpu, '-e', '0.01',
                 '-r', self.pathfinder['dir.cla.scop.txt'], '-m', self.pathfinder['model.tab'],
                 '-p', self.pathfinder['pdbj95d'], '-s', self.pathfinder['self_hits.tab'], self.fastafile,
                 self.tmpfname + '.supfam_scan', self.tmpfname + '.supfam_ass'])

        protein2gdomain = {}
        with open(self.tmpfname + '.supfam_ass') as fin:
            for line in fin:
                accession, match_regions, score, sf_id = \
                    operator.itemgetter(0, 2, 6, 8)([x.strip() for x in line.strip().split('\t')])
                if sf_id == '52592':  # G protein family
                    try:
                        float(score)  # can happen if the entry is '-' (see Echinops telfairi),
                        # in this case, do not consider a valid domain
                        protein2gdomain[accession] = match_regions.split(',')
                    except ValueError:
                        continue
        return protein2gdomain

    def find_best_hit(self):
        # make sure that sequence hits a Rab below specified threshold
        run_cmd([self.pathfinder['phmmer'], '-E', str(self.bh_evalue), '--tblout', self.tmpfname + '.rabs',
                 '--cpu', self.cpu, self.fastafile, self.rabs_db])
        # make sure that sequence hits a Rab at a lower e-value than a non-Rab
        run_cmd([self.pathfinder['phmmer'], '-E', str(self.bh_evalue), '--tblout', self.tmpfname + '.non_rabs',
                 '--cpu', self.cpu, self.fastafile, self.non_rabs_db])

        protein2score = {k: {'evalue_bh_rabs': 10, 'evalue_bh_non_rabs': 10} for k in self.query_id_list}
        with open(self.tmpfname + '.rabs') as fin:
            for line in fin:
                if not line.startswith('#'):
                    query, evalue = operator.itemgetter(2, 4)(line.strip().split())
                    protein2score[query]['evalue_bh_rabs'] = min(protein2score[query]['evalue_bh_rabs'], float(evalue))
        with open(self.tmpfname + '.non_rabs') as fin:
            for line in fin:
                if not line.startswith('#'):
                    query, evalue = operator.itemgetter(2, 4)(line.strip().split())
                    protein2score[query]['evalue_bh_non_rabs'] = min(protein2score[query]['evalue_bh_non_rabs'], float(evalue))
        # If no hits for the threshold set evalue to None
        for v in protein2score.values():
            if v['evalue_bh_rabs'] == 10:
                v['evalue_bh_rabs'] = None
            if v['evalue_bh_non_rabs'] == 10:
                v['evalue_bh_non_rabs'] = None
        return protein2score

    def find_rabf_motifs(self):
        run_cmd([self.pathfinder['mast'], self.rabf_motifs, self.fastafile, '-mt', str(self.motif_evalue), '-hit_list'],
                out=self.tmpfname + '.mast')
        protein2motifs = {k: [] for k in self.query_id_list}
        with open(self.tmpfname + '.mast') as fin:
            for line in fin:
                if not line.startswith('#'):
                    query, motifn, start, end, pvalue = operator.itemgetter(0, 1, 2, 3, 5)(line.strip().split())
                    protein2motifs[query].append([int(motifn), int(start), int(end), float(pvalue)])
        return protein2motifs

    def summarize(self):
        for prot_id, putative_rab in self.gproteins.items():
            print(prot_id, putative_rab.is_rab())

    def _update_fasta(self, ids):
        n = 0
        with open(self.tmpfname, 'w') as fout:
            for seq in SeqIO.parse(self.fastafile, 'fasta'):
                if seq.id in ids:
                    fout.write(seq.format('fasta'))
                    n += 1
        self.fastafile = self.tmpfname
        return n

    def write(self):
        """Write sequences predicted to be Rabs as a fasta file.

        :return: Number of written sequences
        :rtype: int
        """

        rabs = [x.seqrecord for x in self.gproteins.values() if x.is_rab()]
        return SeqIO.write(rabs, self.tmpfname + '.phase2', 'fasta')


class Phase2(object):
    """Assigns Rab subfamilies given Rab protein sequences.

    """

    def __init__(self, fastafile, tmp, **kwargs):
        self.fastafile = fastafile
        self.check_fastafile()

        # Seed database
        def get_db_file(name):
            path = os.path.join(kwargs.get('seed_path', config['seed']['path']), config['seed'][name])
            if os.path.exists(path):
                return path
            else:
                raise IOError("File '{}' doesn't exist".format(path))
        self.rabs_db = kwargs.get('rab_db', get_db_file('rab_db'))
        self.rab_models = kwargs.get('rab_subfamily_hmm', get_db_file('rab_subfamily_hmm'))
        for f in (self.rabs_db, self.rab_models):
            if not os.path.exists(f):
                logger.error('Database file {} is missing.'.format(f))
                raise RuntimeError('Database file {} is missing.'.format(f))
        with open(get_db_file('rab_subfamily_model')) as fin:  # load necessary distributions
            self.subfam2logistic = json.load(fin)

        # Technical parameters
        self.cpu = kwargs.get('cpu', config['param']['cpu'])

        # Classification parameters
        self.bh_evalue = kwargs.get('bh_evalue', config['param']['bh_evalue'])
        self.subfamily_identity = kwargs.get('subfamily_identity', config['param']['subfamily_identity']) * 100
        self.subfamily_score = kwargs.get('subfamily_score', config['param']['subfamily_score'])

        self.pathfinder = Pathfinder(True)
        # Paths for temporary files
        self.tmp = tmp
        self.tmpfname = os.path.join(tmp, os.path.basename(fastafile))  # base name for output files

    def check_fastafile(self):
        if os.stat(self.fastafile).st_size == 0:
            raise IOError('File {} is empty.'.format(self.fastafile))

    def __call__(self):
        """Combine the information generated per subfamily into a single score.

        """

        protein2subfam = {}

        for seq in SeqIO.parse(self.fastafile, 'fasta'):
            protein2subfam[seq.id] = {}
            for rab_subfam in self.subfam2logistic.keys():
                protein2subfam[seq.id][rab_subfam] = {'evalue_hmmscan': 10, 'evalue_phmmer': 10,
                                                      'bh_phmmer': None, 'identity': 0.0}

        #  Run hmmscan on every input protein in .phase2 file
        #  Given an input sequence assign the lowest evalue for each RAB subfamily
        run_cmd([self.pathfinder['hmmscan'], '--tblout', self.tmpfname + '.hmmscan', '--cpu', str(self.cpu),
                 self.rab_models, self.fastafile])
        with open(self.tmpfname + '.hmmscan') as fin:
            for line in fin:
                if not line.startswith('#'):
                    target, query, evalue = operator.itemgetter(0, 2, 4)(line.strip().split())
                    evalue = float(evalue)
                    if evalue < protein2subfam[query][target]['evalue_hmmscan']:
                        protein2subfam[query][target]['evalue_hmmscan'] = evalue

        # Get best hit evalue within each subfamily, try to use phmmer output from Phase1 to speedup calculations
        phmmer_rabs_file = os.path.splitext(self.tmpfname)[0] + '.rabs'
        run_cmd_if_file_missing([self.pathfinder['phmmer'], '-E', str(self.bh_evalue), '--tblout', phmmer_rabs_file,
                                 '--cpu', str(self.cpu), self.fastafile, self.rabs_db], phmmer_rabs_file)
        with open(phmmer_rabs_file) as fin:
            for line in fin:
                if not line.startswith('#'):
                    target, query, evalue_sf = operator.itemgetter(0, 2, 4)(line.strip().split())
                    evalue_sf = float(evalue_sf)
                    sf_rab_subfamily = target.split('___')[2]
                    if query in protein2subfam and evalue_sf < protein2subfam[query][sf_rab_subfamily]['evalue_phmmer']:
                        protein2subfam[query][sf_rab_subfamily]['evalue_phmmer'] = evalue_sf
                        protein2subfam[query][sf_rab_subfamily]['bh_phmmer'] = target

        # Compute sequence identity (not provided by phmmer) of the best (lowest evalue) HSP
        rabs_db_seqs = SeqIO.to_dict(SeqIO.parse(self.rabs_db, 'fasta'))
        for seq in SeqIO.parse(self.fastafile, 'fasta'):
            rs_seq_all = []
            for rab_subfamily in (rs for rs in protein2subfam[seq.id] if protein2subfam[seq.id][rs]['bh_phmmer']):
                rs_seq = rabs_db_seqs[protein2subfam[seq.id][rab_subfamily]['bh_phmmer']]
                rs_seq.id = rab_subfamily
                rs_seq_all.append(rs_seq)
            tmpf_target = tempfile.NamedTemporaryFile()
            SeqIO.write(rs_seq_all, tmpf_target.name, 'fasta')
            tmpf_target.seek(0)
            with open(os.path.devnull, 'w') as devnull:
                proc = subprocess.Popen([self.pathfinder['blastp'], '-subject', tmpf_target.name, '-max_hsps', '1',
                                         '-outfmt', '10 sseqid pident'],
                                        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=devnull)
                stdout_value = proc.communicate(seq.format('fasta').encode())[0].decode()
            for rs, pident in (line.split(',') for line in stdout_value.strip().split('\n') if line):
                protein2subfam[seq.id][rs]['identity'] = int(float(pident))

        self.score = {}
        self.annotation = {}
        for prot_id in protein2subfam:
            numerators = {}
            denominator = 0.0
            for rab_subfamily in protein2subfam[prot_id]:
                # transform evalues
                try:
                    t_evalue_phmmer = (-1) * math.log10(protein2subfam[prot_id][rab_subfamily]['evalue_phmmer'])
                except ValueError:
                    t_evalue_phmmer = 350
                try:
                    t_evalue_hmmscan = (-1) * math.log10(protein2subfam[prot_id][rab_subfamily]['evalue_hmmscan'])
                except ValueError:
                    t_evalue_hmmscan = 350
                # transform log(p-value) into value [0, 1] according to sigmoidal curve
                x1 = scipy.stats.logistic.cdf(t_evalue_phmmer, self.subfam2logistic[rab_subfamily]['ph_location'],
                                              self.subfam2logistic[rab_subfamily]['ph_scale'])
                x2 = scipy.stats.logistic.cdf(t_evalue_hmmscan, self.subfam2logistic[rab_subfamily]['hs_location'],
                                              self.subfam2logistic[rab_subfamily]['ph_scale'])
                # complete naive bayesian classifier computations; 1e+12 is just to not run into numerical
                # problems because of too small values, cancels out after division
                numerators[rab_subfamily] = 1e+12 * x1 * x2
                denominator += numerators[rab_subfamily]
            #compute the final score
            for rab_subfamily in protein2subfam[prot_id]:
                self.score.setdefault(prot_id, {})[rab_subfamily] = numerators[rab_subfamily] / denominator

            # get values to test the remaining two conditions that may lead to RabXs
            # FIXME is this correct, when evalue_phmmer is identical for few subfamilies the one with lowest identity
            # will be taken
            seq_identity_of_best_hit = min((rab_subfamily['evalue_phmmer'], rab_subfamily['identity'])
                                           for rab_subfamily in protein2subfam[prot_id].values())[1]
            max_subfamily_score, rab_subfamily_predicted = max((rab_subfamily_score, rab_subfamily)
                                                               for rab_subfamily, rab_subfamily_score
                                                               in self.score[prot_id].items())
            if seq_identity_of_best_hit <= self.subfamily_identity or max_subfamily_score < self.subfamily_score:
                for rab_subfamily in self.score[prot_id]:
                    self.score[prot_id][rab_subfamily] /= 2
                rab_subfamily_predicted = 'rabX'
                self.score[prot_id]['rabX'] = 0.5

            self.annotation[prot_id] = (rab_subfamily_predicted, self.score[prot_id][rab_subfamily_predicted])

    def get_best_subfamily(self):
        pass

    def get_top_subfamilies(self):
        pass

    def summarize(self):
        print(self.annotation)


class Rabmyfire(object):
    """ Classifies proteins as (non)Rabs, assigns subfamily if applicable.
    """

    def __init__(self, **kwargs):
        self.params = kwargs
        self.tmp = tempfile.mkdtemp(prefix='rabifier_', dir=config['tmp'])
        self.check()

    def check(self):
        """ Check if data and third party tools, necessary to run the classification, are available

        :raises: RuntimeError
        """

        pathfinder = Pathfinder(True)
        if pathfinder.add_path(pathfinder['superfamily']) is None:
            raise RuntimeError("'superfamily' data directory is missing")

        for tool in ('hmmscan', 'phmmer', 'mast', 'blastp', 'ass3.pl', 'hmmscan.pl'):
            if not pathfinder.exists(tool):
                raise RuntimeError('Dependency {} is missing'.format(tool))

    def __call__(self, fastafile):
        """ Run the classifier

        :param str fastafile: path to a fasta file with protein sequences
        :return: rab predictions
        :rtype: iter(Gprotein)
        """

        try:
            phase1 = Phase1(os.path.abspath(fastafile), self.tmp, **self.params)
            returncode, err_msg = phase1()
            if returncode != 0:
                raise RuntimeError(err_msg)
            phase1.write()

            try:
                phase2 = Phase2(phase1.tmpfname + '.phase2', self.tmp, **self.params)
                phase2()

                for protein_name, putative_rab in phase1.gproteins.items():
                    if protein_name in phase2.annotation:  # Update annotation if protein was in Phase2
                        putative_rab.rab_subfamily, putative_rab.rab_subfamily_score = phase2.annotation[protein_name]
                        putative_rab.rab_subfamily_top5 = \
                            sorted(phase2.score[protein_name].items(), key=lambda x: x[1], reverse=True)[:5]
            except IOError:
                pass
        finally:
            shutil.rmtree(self.tmp)
            pass

        return phase1.gproteins.values()


def main():
    parser = argparse.ArgumentParser(description="Rabifier: a bioinformatic classifier of Rab GTPases (v{})".format(__version__))
    parser.add_argument('-v', '--version', action='version', version=__version__)
    parser.add_argument('input', help="input file name, a fasta file containing protein sequence(s)")
    parser.add_argument('-o', '--output', help="output file name [-]", type=argparse.FileType('w'), default='-')
    parser.add_argument('--outfmt', help="output format [text]", choices=('text', 'json', 'csv'), default='text')
    parser.add_argument('--show_positive', action='store_true', help="include only Rab positive predictions in the "
                                                                      "output")
    parser.add_argument('--cpu', help="maximal number of threads to use [{}]".format(config['param']['cpu']), type=int,
                        default=config['param']['cpu'])
    parser.add_argument('--fast', action='store_true', help="phase 1 speedup: stop if it is not a G protein "
                                                            "(does not show the full summary)")
    parser.add_argument('--bh_evalue', type=float, default=config['param']['bh_evalue'],
                        help="e-value threshold for best hit search [{:.0e}]".format(config['param']['bh_evalue']))
    parser.add_argument('--motif_evalue', type=float, default=config['param']['motif_evalue'],
                        help="e-value threshold for motif search [{:.0e}]".format(config['param']['motif_evalue']))
    parser.add_argument('--motif_number', type=int, default=config['param']['motif_number'],
                        help="minimum number of RabF motifs [{}]".format(config['param']['motif_number']))
    parser.add_argument('--subfamily_identity', type=float, default=config['param']['subfamily_identity'],
                        help="minimum sequence identity with subfamily's best hit [{}]".format(
                            config['param']['subfamily_identity']))
    parser.add_argument('--subfamily_score', type=float, default=config['param']['subfamily_score'],
                        help="minimum subfamily score [{}]".format(config['param']['subfamily_score']))

    try:
        args = parser.parse_args()
        logger.debug(args)
        open(args.input).close()  # check if the file exists
        try:
            classifier = Rabmyfire(cpu=args.cpu,
                                   fast=args.fast,
                                   bh_evalue=args.bh_evalue,
                                   motif_evalue=args.motif_evalue,
                                   motif_number=args.motif_number,
                                   subfamily_identity=args.subfamily_identity,
                                   subfamily_score=args.subfamily_score
                                   )
            predictions = classifier(args.input)
            if args.outfmt == 'text':
                for putative_rab in predictions:
                    if putative_rab.is_rab() or not args.show_positive:
                        args.output.write(str(putative_rab) + '\n')
            elif args.outfmt == 'json':
                d = {putative_rab.seqrecord.id: putative_rab.to_dict() for putative_rab in predictions
                     if putative_rab.is_rab() or not args.show_positive}
                json.dump(d, args.output, indent=2)
            elif args.outfmt == 'csv':
                writer = csv.writer(args.output)
                writer.writerow(['id', 'g_domain', 'evalue_bh_rabs', 'evalue_bh_nonrabs', 'rabf', 'is_rab',
                                 'subfamily', 'top_5'])
                writer.writerows([putative_rab.to_list() for putative_rab in predictions
                                  if putative_rab.is_rab() or not args.show_positive])
        except IOError as err:
            parser.error(str(err))
        except RuntimeError as err:
            parser.error(str(err))
        finally:
            args.output.close()

    except IOError as err:
        print(err.args, err.filename, err.errno)
        parser.error(str(err))
