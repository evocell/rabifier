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

from __future__ import print_function

import os
import json
from collections import defaultdict
import math
import operator
import logging
import argparse

from Bio import SeqIO, AlignIO
import numpy as np
import scipy.stats
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

from . import __version__
from . import config
from .utils import Pathfinder, run_cmd, run_cmd_if_file_missing, merge_files


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Seed(object):
    """Seed requires data produced by fetchdb.py

    FIXME temporary files and directories are not removed at exit
    FIXME final files are not copied to the destination
    """

    def __init__(self, source, **kwargs):
        self.source = os.path.abspath(source)
        self.pathfinder = Pathfinder(True)

        # Temporary folders, create if do not exist
        self.path_tmp = kwargs.get('tmp', os.path.join(config['tmp'], 'rabifier_tmp'))
        for name in ('build', 'test', 'seed'):
            path = os.path.join(self.path_tmp, name)
            if not os.path.exists(path):
                os.makedirs(path)

        self.path = {
            'non_rab_db': os.path.join(self.source, 'other', '{}.full'.format(config['seed']['non_rab_db'])),
            'rab_subfamily_model_manual_override':
                os.path.join(self.source, 'other', 'rab_subfamilies_logreg_params.dat'),
            'rab_subfamily': os.path.join(self.source, 'rab', '{}.fasta.full')
        }

        self.output = {
            'rab_db': os.path.join(self.path_tmp, 'seed', config['seed']['rab_db']),
            'rab_subfamily_hmm': os.path.join(self.path_tmp, 'seed', config['seed']['rab_subfamily_hmm']),
            'rab_subfamily_model': os.path.join(self.path_tmp, 'seed', config['seed']['rab_subfamily_model']),
            'rab_f': os.path.join(self.path_tmp, 'seed', config['seed']['rab_f']),
            'non_rab_db': os.path.join(self.path_tmp, 'seed', config['seed']['non_rab_db']),
            'rab_db_reduced': os.path.join(self.path_tmp, 'build', 'rab_db_reduced.fasta'),
            'rab_db_reduced_msa': os.path.join(self.path_tmp, 'build', 'rab_db_reduced.mafft'),
            'rab_f_single_motif': os.path.join(self.path_tmp, 'build', 'motif_{}.meme')
        }

        # Technical parameters
        self.cpu = str(kwargs.get('cpu', config['param']['cpu']))

        # Check if the datasets and third party tools, necessary to prepare the seed database, are available
        self.check()

        # Load subfamily names that will be considered for the seed database construction
        self.rab_subfamilies = config['rab_subfamilies']

    def get_subfamily_path(self, subfamily, extension=None):
        path = os.path.join(self.path_tmp, 'build', subfamily)
        if extension is None:
            return path
        else:
            return '.'.join([path, extension])

    def check(self):
        """ Check if data and third party tools are available

        :raises: RuntimeError
        """

        #for path in self.path.values():
        #    if not os.path.exists(path):
        #        raise RuntimeError("File '{}' is missing".format(path))

        for tool in ('cd-hit', 'prank', 'hmmbuild', 'hmmpress', 'hmmscan', 'phmmer', 'mafft', 'meme'):
            if not self.pathfinder.exists(tool):
                raise RuntimeError("Dependency {} is missing".format(tool))

    def __call__(self):
        # NOTE For some long running commands skip execution if the output already exists. Not sure if it completely safe.
        self.build_subfamilies()
        self.build_subfamily_models()
        self.build_rab_f_models()
        self.generate_non_rabs()

    def dump_db_files(self, destination):
        """Copy final (i.e. necessary to annotate new Rabs) database files to the destination.

        :param str destination: Path to the destination directory
        :raise IOError: Raises IOError if any database file does not exist
        """
        pass

    def build_subfamilies(self):
        for subfamily in self.rab_subfamilies:
            logging.info('Generating {} subfamily pHMM'.format(subfamily))
            run_cmd([self.pathfinder['cd-hit'], '-i', self.path['rab_subfamily'].format(subfamily),
                     '-o', self.get_subfamily_path(subfamily, 'fasta'), '-d', '100',
                     '-c', str(config['param']['subfamily_identity_threshold']), '-g', '1', '-T', self.cpu])

            if not os.path.exists(self.get_subfamily_path(subfamily, 'sto')):
                if sum(1 for _ in SeqIO.parse(self.get_subfamily_path(subfamily, 'fasta'), 'fasta')) > 1:
                    # NOTE here we try to reuse old alignment to skip long computation
                    run_cmd_if_file_missing([self.pathfinder['prank'],
                                             '-d={}'.format(self.get_subfamily_path(subfamily, 'fasta')),
                                             '-o={}'.format(self.get_subfamily_path(subfamily)),
                                             '-F'],
                                            self.get_subfamily_path(subfamily, 'sto'))
                    AlignIO.convert(self.get_subfamily_path(subfamily, 'best.fas'), 'fasta',
                                    self.get_subfamily_path(subfamily, 'sto'), 'stockholm')
                else:
                    # No MSA if only one sequence
                    SeqIO.convert(self.get_subfamily_path(subfamily, 'fasta'), 'fasta',
                                  self.get_subfamily_path(subfamily, 'sto'), 'stockholm')

            run_cmd([self.pathfinder['hmmbuild'], '-n', subfamily,
                     '-O', self.get_subfamily_path(subfamily, 'sto.hmmbuild'),
                     '--amino', '--cpu', self.cpu, self.get_subfamily_path(subfamily, 'hmm'),
                     self.get_subfamily_path(subfamily, 'sto')])

        logging.info('Building Rab subfamilies pHMM DB')

        merge_files([self.get_subfamily_path(subfamily, 'hmm') for subfamily in self.rab_subfamilies],
                    self.output['rab_subfamily_hmm'])
        run_cmd([self.pathfinder['hmmpress'], '-f', self.output['rab_subfamily_hmm']])

        merge_files([self.get_subfamily_path(subfamily, 'fasta') for subfamily in self.rab_subfamilies],
                    self.output['rab_db'])

    def build_subfamily_models(self):

        param_subfam_score_top_threshold = 10
        param_subfam_score_ph_location_offset = 30
        param_subfam_score_ph_scale = 3
        param_subfam_score_hs_location_offset = 40
        param_subfam_score_hs_scale = 5
        param_subfam_score_hs_location_offset_small = 12

        def sort_replace_zeros_truncate(values, truncate):
            values.sort()
            # defines the default value if all values are 0.0 !!!
            # FIXME this seems to work as the minimal float, yet not sure if works for different
            # Python versions/machines.
            minimal = 1e-323
            for x in values:
                if x > minimal:
                    minimal = x
                    break
            if truncate:
                return [x if x > 0 else minimal for x in values][:param_subfam_score_top_threshold]
            else:
                return [x if x > 0 else minimal for x in values]

        rab_subfam2evalues_phmmer = {}
        rab_subfam2evalues_phmmer_self = defaultdict(list)
        rab_subfam2evalues_hmmscan = defaultdict(list)

        for subfamily in self.rab_subfamilies:
            logging.info('Computing E values for {} subfamily'.format(subfamily))
            # Run phmmer for each sequence in a subfamily against database of sequences from all subfamilies
            # NOTE '-E', E_VALUE_THRESHOLD should be omitted (If not, good to spot misannotations)
            run_cmd_if_file_missing([self.pathfinder['phmmer'],
                                     '--tblout', self.get_subfamily_path(subfamily, 'phmmer'),
                                     '--cpu', self.cpu,
                                     self.get_subfamily_path(subfamily, 'fasta'),
                                     self.output['rab_db']],
                                    self.get_subfamily_path(subfamily, 'phmmer'))

            query2evalues = defaultdict(list)
            with open(self.get_subfamily_path(subfamily, 'phmmer')) as fin:
                for line in fin:
                    if not line.startswith('#'):
                        target, query, evalue = operator.itemgetter(0, 2, 4)(line.strip().split())
                        # DO NOT COUNT THE ALIGNMENT TO ITSELF !!!!! CHECK THE CONSEQUENCES OF THAT !!!
                        if target.split('___')[2] == subfamily:
                            if target == query:
                                rab_subfam2evalues_phmmer_self[subfamily].append(float(evalue))
                            else:
                                query2evalues[query].append(float(evalue))

            rab_subfam2evalues_phmmer[subfamily] = []
            for query, evalues in query2evalues.items():
                # for each subfamily one has # members(==#queries) * param_subfam_score_top_threshold evalues!
                rab_subfam2evalues_phmmer[subfamily].extend(sort_replace_zeros_truncate(evalues, True))

            # Run hmmscan for each sequence in a subfamily against a rab subfamily profile database
            run_cmd_if_file_missing([self.pathfinder['hmmscan'],
                                     '--tblout', self.get_subfamily_path(subfamily, 'hmmscan'),
                                     '--cpu', self.cpu,
                                     self.output['rab_subfamily_hmm'],
                                     self.get_subfamily_path(subfamily, 'fasta')],
                                    self.get_subfamily_path(subfamily, 'hmmscan'))

            with open(self.get_subfamily_path(subfamily, 'hmmscan')) as fin:
                for line in fin:
                    if not line.startswith('#'):
                        target, query, evalue = operator.itemgetter(0, 2, 4)(line.strip().split())
                        if target == subfamily:
                            rab_subfam2evalues_hmmscan[subfamily].append(float(evalue))
            rab_subfam2evalues_hmmscan[subfamily] = \
                sort_replace_zeros_truncate(rab_subfam2evalues_hmmscan[subfamily], truncate=False)

        logging.info('Fitting distributions to subfamily E values')

        rab_subfamily2logistic_fit = {}

        # load manual logistic fit parameters, which override parameters estimated from the data
        # TODO I tried to fit one per subfamily, but then tings like favor Rab31, ypt11, i.e. highly
        # skewed distributions !!!! ONLY IF VALUES ARE CONSISTENTLY LOW !!!!! hope to catch those cases by
        # <40%ID etc. !!
        # EXCEPTION !! otherwise fit is too sharp !!
        manual_params = {}
        with open(self.path['rab_subfamily_model_manual_override']) as fin:
            for line in fin:
                l = line.strip().split()
                if l and not l[0].startswith('#'):
                    manual_params[l[0]] = [float(x) if x[0].isdigit() else None for x in l[1:]]

        for subfamily in self.rab_subfamilies:
            ev_ph = np.log10(rab_subfam2evalues_phmmer[subfamily]) * -1
            if subfamily in manual_params and all(manual_params[subfamily][0:2]):
                ph_location, ph_scale = manual_params[subfamily][0:2]
            else:
                if len(ev_ph) > 1:
                    ph_location, ph_scale = scipy.stats.logistic.fit(ev_ph)
                else:
                    # Override ev_ph
                    ev_ph = np.log10(rab_subfam2evalues_phmmer_self[subfamily]) * -1
                    ph_location = ev_ph[0] - param_subfam_score_ph_location_offset
                    ph_scale = param_subfam_score_ph_scale
            rab_subfamily2logistic_fit[subfamily] = {'ph_location': ph_location, 'ph_scale': ph_scale}

            ev_hs = np.log10(rab_subfam2evalues_hmmscan[subfamily]) * -1

            if subfamily in manual_params and all(manual_params[subfamily][2:4]):
                hs_location, hs_scale = manual_params[subfamily][2:4]
            else:
                if len(ev_hs) > 1:
                    hs_location, hs_scale = scipy.stats.logistic.fit(ev_hs)
                    # just to make sure that the majority of the values lies beneath the curve !!!
                    hs_location -= param_subfam_score_hs_location_offset_small
                else:
                    hs_location = ev_hs[0] - param_subfam_score_hs_location_offset
                    hs_scale = param_subfam_score_hs_scale
            rab_subfamily2logistic_fit[subfamily]['hs_location'] = hs_location
            rab_subfamily2logistic_fit[subfamily]['hs_scale'] = hs_scale

            if plt is not None:
                fig, axarr = plt.subplots(2, 2)
                fig.subplots_adjust(hspace=0.3)
                fig.suptitle(subfamily)
                for ax, data, method in ((axarr[0, 0], ev_ph, 'phmmer'), (axarr[1, 0], ev_hs, 'hmmscan')):
                    ax.hist(data, bins=20, range=(min(data) - 10, max(data) + 10), color='#eeeeee', histtype='bar',
                            rwidth=0.7)
                    ax.set_xlabel(method + r' $(-1) \log_{10}$(e-value)')
                    ax.set_ylabel('count')

                for ax, data, location, scale, method in ((axarr[0, 1], ev_ph, ph_location, ph_scale, 'phmmer'),
                                                          (axarr[1, 1], ev_hs, hs_location, hs_scale, 'hmmscan')):
                    ax.hist(data, bins=50, range=(0, 325), color='#eeeeee', cumulative=True, histtype='stepfilled')
                    x, y = zip(*[(i, scipy.stats.logistic.cdf(i, location, scale) * len(data)) for i in range(325)])
                    ax.plot(x, y, lw=2, color='red')
                    ax.text(10, 0.8 * len(data),
                            'logistic fit\nloc: {:.1f}\nscale: {:.1f}'.format(location, scale), fontsize=9)
                    ax.set_xlabel(method + r' $(-1) \log_{10}$(e-value)')

                plt.savefig(self.get_subfamily_path(subfamily + '_model', 'pdf'))
                plt.close()

        with open(self.output['rab_subfamily_model'], 'w') as fout:
            json.dump(rab_subfamily2logistic_fit, fout, indent=2)

    def test_reference_set(self):
        """ Test for outliers in the reference data set.
        """

        logging.info('Testing for outliers in the reference set')
        with open(os.path.join(self.path_tmp, 'test', 'summary.txt'), 'w') as fout:
            for subfamily in self.rab_subfamilies:
                fout.write('=== {} ===\n\n'.format(subfamily))

                testfile = os.path.join(self.path_tmp, 'test', subfamily + '.phmmer')
                run_cmd_if_file_missing([self.pathfinder['phmmer'], '--tblout', testfile, '--cpu', self.cpu,
                                         self.get_subfamily_path(subfamily, 'fasta'),
                                         self.get_subfamily_path(subfamily, 'fasta')],
                                        testfile)
                query2evalues = defaultdict(list)
                evalue_query_target = []
                with open(testfile) as fin:
                    for line in fin:
                        if not line.startswith('#'):
                            target, query, evalue = operator.itemgetter(0, 2, 4)(line.strip().split())
                            try:
                                evalue = -1 * math.log10(float(evalue))
                            except ValueError:
                                evalue = 350.0
                            evalue_query_target.append((evalue, query, target))
                            query2evalues[query].append(evalue)

                for query, evalues in query2evalues.items():
                    query2evalues[query] = np.mean(evalues)

                fout.write('# Worst average scores\n')
                for a, b in sorted(query2evalues.items(), key=operator.itemgetter(1))[:15]:
                    fout.write('{:.2f}\t{}\n'.format(b, a))
                fout.write('\n')

                fout.write('# Worst pairwise scores\n')
                for a, b, c in sorted(evalue_query_target)[:15]:
                    fout.write('{:.2f}\t{:25s}\t{}\n'.format(a, b, c))
                fout.write('\n')

                run_cmd_if_file_missing([self.pathfinder['hmmpress'], self.get_subfamily_path(subfamily, 'hmm')],
                                        self.get_subfamily_path(subfamily, 'hmm.h3m'))
                testfile = os.path.join(self.path_tmp, 'test', subfamily + '.hmmscan')
                run_cmd_if_file_missing([self.pathfinder['hmmscan'], '--tblout', testfile, '--cpu', self.cpu,
                                         self.get_subfamily_path(subfamily, 'hmm'),
                                         self.get_subfamily_path(subfamily, 'fasta')],
                                        testfile)
                evalue_query = []
                with open(testfile) as fin:
                    for line in fin:
                        if not line.startswith('#'):
                            query, evalue = operator.itemgetter(2, 4)(line.strip().split())
                            try:
                                evalue_query.append(((-1)*math.log10(float(evalue)), query))
                            except ValueError:
                                evalue_query.append((350.0, query))

                fout.write('# Worst scores against the subfamily HMM\n')
                for a, b in sorted(evalue_query)[:15]:
                    fout.write('{:.2f}\t{}\n'.format(a, b))
                fout.write('\n')

    def build_rab_f_models(self):
        motifs = ('IGVDF', 'KLQIW', 'RFxxxT', 'YYRGA', 'LVYDIT')

        logging.info('Building RabF motif models')
        # Reduce redundancy
        run_cmd_if_file_missing([self.pathfinder['cd-hit'], '-i', self.output['rab_db'], '-o',
                                 self.output['rab_db_reduced'], '-d', '100',
                                 '-c', str(config['param']['rab_f_identity_threshold']), '-g', '1'],
                                self.output['rab_db_reduced'])
        run_cmd_if_file_missing([self.pathfinder['mafft'], '--thread', self.cpu, self.output['rab_db_reduced']],
                                self.output['rab_db_reduced_msa'], out=self.output['rab_db_reduced_msa'])

        # Remove gapped columns from the MSA
        alignment = AlignIO.read(self.output['rab_db_reduced_msa'], 'fasta')
        keep = []
        for i in range(alignment.get_alignment_length()):
            symbols = set(alignment[:, i]) - {'-', 'X', 'x'}
            # Keep column if at least two sequences have informative residues at this position
            if len(symbols) > 1 or (len(symbols) == 1 and alignment[:, i].count(symbols.pop()) > 1):
                keep.append(i)
        trimmed_alignment = alignment[:, keep[0]:keep[0]+1]
        for i in keep[1:]:
            trimmed_alignment += alignment[:, i:i+1]
        AlignIO.write([trimmed_alignment], self.output['rab_db_reduced_msa'] + '.trimmed', 'fasta')

        # Recast MSA to unaligned sequences
        records = []
        for record in SeqIO.parse(self.output['rab_db_reduced_msa'] + '.trimmed', 'fasta'):
            record.seq = record.seq.ungap(gap='-')
            records.append(record)
        SeqIO.write(records, self.output['rab_db_reduced_msa'] + '.trimmed.clean', 'fasta')

        # Build models for all motifs
        for motif in motifs:
            run_cmd_if_file_missing([self.pathfinder['meme'], self.output['rab_db_reduced_msa'] + '.trimmed.clean',
                                     '-protein', '-cons', motif, '-w', '{:d}'.format(len(motif)), '-mod', 'zoops',
                                     '-text', '-maxsize', '1000000', '-p', self.cpu],
                                    self.output['rab_f_single_motif'].format(motif),
                                    out=self.output['rab_f_single_motif'].format(motif))

        # Parse into MEME Minimal Motif Format
        with open(self.output['rab_f'], 'w') as fout:
            # Header
            fout.write('MEME version 4\n\nALPHABET= ACDEFGHIKLMNPQRSTVWY\n\n')
            # Background letter frequencies
            background_freqs = None
            with open(self.output['rab_f_single_motif'].format(motifs[0])) as fin:
                for line in fin:
                    if line.startswith('Background letter frequencies'):
                        background_freqs = next(fin) + next(fin) + next(fin)
                        break

            fout.write('Background letter frequencies\n{}\n'.format(background_freqs))
            # Motifs
            for motif in motifs:
                fout.write('MOTIF {}\n'.format(motif))
                with open(self.output['rab_f_single_motif'].format(motif)) as fin:
                    for line in fin:
                        if line.startswith('letter-probability matrix: '):
                            fout.write(line)
                            for i in range(len(motif)):
                                fout.write(next(fin).lstrip())
                            fout.write('\n')
                            break

    def generate_non_rabs(self):
        """ Shrink the non-Rab DB size by reducing sequence redundancy.
        """

        logging.info('Building non-Rab DB')
        run_cmd([self.pathfinder['cd-hit'], '-i', self.path['non_rab_db'], '-o', self.output['non_rab_db'],
                 '-d', '100', '-c', str(config['param']['non_rab_db_identity_threshold']), '-g', '1', '-T', self.cpu])
        os.remove(self.output['non_rab_db'] + '.clstr')


def main():
    parser = argparse.ArgumentParser(description="Rabifier: a bioinformatic classifier of Rab GTPases (v{})".format(__version__))
    parser.add_argument('-v', '--version', action='version', version=__version__)
    parser.add_argument('source', help="Source data directory", type=str)
    parser.add_argument('-o', '--output', help="path to a directory where the output files will be stored [{}]".format(
        os.path.join(config['tmp'], 'rabifier')), type=str, default=os.path.join(config['tmp'], 'rabifier'))
    parser.add_argument('--cpu', help="maximal number of threads to use [{}]".format(config['param']['cpu']), type=int,
                        default=config['param']['cpu'])

    args = parser.parse_args()
    seed = Seed(source=args.source, cpu=args.cpu, tmp=args.output)
    seed()
    seed.test_reference_set()
