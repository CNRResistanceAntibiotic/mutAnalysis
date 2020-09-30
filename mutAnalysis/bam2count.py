#!/usr/bin/env python3
"""
The version is stored here in a separate file so it can exist in only one place.
https://stackoverflow.com/a/7071358/2438989
Copyright 2020 Aurélien BIRER (abirer36@gmail.com)
https://github.com/CNRResistanceAntibiotic/mutAnalysis

This file is part of CGST. CGST is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. CGST is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with CGST. If
not, see <http://www.gnu.org/licenses/>.
"""
import os
import re
from collections import OrderedDict
from subprocess import Popen, PIPE, STDOUT
import numpy as np
import pandas as pd


def write_site_file(position, output_dir):
    pattern = re.compile('([0-9a-zA-Z_-]+):*([0-9]*)-*([0-9]*)')
    match = pattern.match(position)
    if match:
        txt = ''
        site_file = os.path.join(output_dir, 'site_file.txt')
        for item in match.groups():
            if item != '':
                txt = txt + '\t{0}'.format(item)
        with open(site_file, 'w') as out_f:
            out_f.write(txt[1:])
        return site_file
    else:
        print('\nPositions not identified\n')
        print('Position formats:\n'
              '\t<contig or chromosome name>\n'
              '\t<contig or chromosome name>:<position in chromome or contig>\n'
              '\t<contig or chromosome name>:<start position in chromome or contig>:<end position in chromome or contig>')
        exit(1)


def bam_count(bam_file, fasta_ref, output_dir, q=0, b=0, feature_name='', site_file='', force=False):
    sample = os.path.basename(bam_file).split(".")[0]
    if feature_name == '' and site_file == '':
        out_file = os.path.join(output_dir, '{0}_raw.csv'.format(sample))
        cmd = '$(which bam-readcount) -w 0 -q {0} -b {1} -i -f {2} {3} > {4}'.format(q, b, fasta_ref, bam_file,
                                                                                     out_file)
    elif feature_name != '' and site_file == '':
        out_file = os.path.join(output_dir, '{0}_{1}_raw.csv'.format(sample, feature_name))
        cmd = '$(which bam-readcount) -w 0 -q {0} -b {1} -i -f {2} {3} > {4}'.format(q, b, fasta_ref, bam_file,
                                                                                     out_file)
    else:
        out_file = os.path.join(output_dir, '{0}_{1}_raw.csv'.format(sample, feature_name))
        cmd = '$(which bam-readcount) -w 0 -q {0} -b {1} -i -l {2} -f {3} {4} > {5}'.format(q, b, site_file, fasta_ref,
                                                                                      bam_file, out_file)
    if not os.path.exists(out_file) or force:
        process = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT).stdout.read()
        log_info = "Command line executed: {0}\n\n\n{1}".format(cmd, process.decode("utf-8"))
        print(log_info)
    else:
        print('\nRead count file {0} already done.\n'.format(out_file))
    return out_file


def calculate_stats(data):
    result = OrderedDict()
    result['perc>=30'] = round(100 * len([x for x in data if x >= 30]) / float(len(data)), 2)
    result['perc>=20'] = round(100 * len([x for x in data if x >= 20]) / float(len(data)), 2)
    result['perc>=10'] = round(100 * len([x for x in data if x >= 10]) / float(len(data)), 2)
    result['mean'] = round(sum(data) / float(len(data)), 2)
    result['50_perc'] = round(np.percentile(np.array(data), 50), 2)
    result['25_perc'] = round(np.percentile(np.array(data), 25), 2)
    result['75_perc'] = round(np.percentile(np.array(data), 75), 2)
    result['max'] = max(data)
    result['min'] = min(data)
    return result


def bam_count_stats(bam_count_file, feature_name, header, output_dir, bam_file):

    sample = os.path.basename(bam_file).split(".")[0]

    if feature_name == '':
        out_file = os.path.join(output_dir, '{0}_{1}_stats'.format(sample, 'whole_genome'))
    else:
        out_file = os.path.join(output_dir, '{0}_{1}_stats'.format(sample, feature_name))
    with open(bam_count_file) as count_f:
        ctgs = []
        result_stat = []
        result_data = []
        start = 0
        end = 0
        base_nb = 0

        for line in count_f:
            line = line.strip().split('\t')
            ctg = line[0]
            if ctg not in ctgs:
                if len(ctgs) > 0:
                    d = OrderedDict([('ID', ctgs[-1]), ('start', start), ('end', end), ('size', base_nb)])
                    result_data.append(d)
                    s = []
                    # for key, value in calculate_stats(ctg_depth).iteritems():
                    #    s.append(('All_depth_{0}'.format(key), value))
                    for key, value in calculate_stats(ctg_ref_depth).items():
                        s.append(('Ref_depth_{0}'.format(key), value))
                    for key, value in calculate_stats(ctg_ref_qual).items():
                        s.append(('Ref_quali_{0}'.format(key), value))
                    s = OrderedDict(s)
                    result_stat.append(s)
                ctg_ref_depth, ctg_ref_qual = [], []
                # ctg_ref_depth = []
                ctgs.append(ctg)
                start = -1
                end = 1
                base_nb = 0

            pos = int(line[1])
            if start == -1:
                start = pos
            if pos >= end:
                end = pos
            base_nb += 1

            ref = line[2]
            # ctg_depth.append(int(line[3]))

            for item in line[4:]:
                data = dict(zip(header.split(':'), item.split(':')))
                base = data['base']
                if base == ref:
                    ctg_ref_depth.append(int(data['count']))
                    ctg_ref_qual.append(round(float(data['avg_base_quality']), 2))

        d = OrderedDict([('ID', ctgs[-1]), ('start', start), ('end', end), ('size', base_nb)])
        result_data.append(d)
        s = []
        # for key, value in calculate_stats(ctg_depth).iteritems():
        #    s.append(('All_depth_{0}'.format(key), value))
        for key, value in calculate_stats(ctg_ref_depth).items():
            s.append(('Ref_depth_{0}'.format(key), value))
        for key, value in calculate_stats(ctg_ref_qual).items():
            s.append(('Ref_quali_{0}'.format(key), value))
        s = OrderedDict(s)
        result_stat.append(s)

        df_data = pd.DataFrame.from_dict(result_data)
        df_stat = pd.DataFrame.from_dict(result_stat)
        weighted_mean = ((df_stat.multiply(df_data['size'], axis=0).sum()).div(df_data['size'].sum(), axis=0)).round(2)
        df_stat = df_stat.append(weighted_mean, ignore_index=True)
        df_stat = df_stat.round(2)
        df_data = df_data.append([{'ID': 'Overall', 'size': df_data['size'].sum(), 'end': '-', 'start': '-'}],
                                 ignore_index=True)
        df = pd.concat([df_data, df_stat], axis=1)
        df.sort_values('size', inplace=True)
        df.to_csv('{0}.csv'.format(out_file), sep='\t', header=True, index=True)
        df.to_html('{0}.html'.format(out_file), index=True)
        return '{0}.csv'.format(out_file)


def bam_count_extract(bam_count_file, feature_name, header, output_dir, bam_file):

    sample = os.path.basename(bam_file).split(".")[0]
    if feature_name == '':
        out_file = os.path.join(output_dir, sample + '_{0}_count'.format('whole_genome'))
    else:
        out_file = os.path.join(output_dir, sample + '_{0}_count'.format(feature_name))

    with open(bam_count_file) as count_f:
        result_data = []
        for line in count_f:
            line = line.strip().split('\t')
            ctg = line[0]
            pos = int(line[1])
            ref = line[2]
            depth = int(line[3])
            data = [('ID', ctg), ('position', pos), ('reference', ref), ('total_depth', depth)]
            for item in line[4:]:
                d = dict(zip(header.split(':'), item.split(':')))
                base = d['base']
                if base != '=':
                    data.append(('{0}_depth'.format(base), int(d['count'])))
                    data.append(('{0}_quality'.format(base), round(float(d['avg_base_quality']), 2)))
            data = OrderedDict(data)
            result_data.append(data)

        df = pd.DataFrame.from_dict(result_data)
        # print(df)
        df.to_csv('{0}.csv'.format(out_file), sep='\t', header=True, index=True)
        df.to_html('{0}.html'.format(out_file), index=True)
    return '{0}.csv'.format(out_file)


def main(wk_dir, sequence_file, position, feature_name="", stats=True, data=True):

    bam_file = os.path.join(wk_dir, 'sequence.bam')

    site_file = ""

    while feature_name == '' and position != '':
        # feature_name = raw_input('Enter a feature name for the position and validate by enter:\n')
        print("YOU NEED TO SELECT A FEATURE")
        exit(1)
    if feature_name == '':
        position = ''
    if position != '':
        site_file = write_site_file(position, wk_dir)

    bam_count_file = bam_count(bam_file, sequence_file, wk_dir, 0, 0, feature_name, site_file, force=False)

    #if os.path.exists(site_file):
    #    os.remove(site_file)

    header = 'base:count:avg_mapping_quality:avg_base_quality:avg_se_mapping_quality:num_plus_strand:num_minus_strand:' \
             'avg_pos_as_fraction:avg_num_mismatches_as_fraction:avg_sum_mismatch_qualities:num_q2_containing_reads:' \
             'avg_distance_to_q2_start_in_q2_reads:avg_clipped_length:avg_distance_to_effective_3p_end'

    out_file_list = [bam_count_file]
    if stats:
        stat_file = bam_count_stats(bam_count_file, feature_name, header, wk_dir, bam_file)
        out_file_list.append(stat_file)
    if data:
        data_file = bam_count_extract(bam_count_file, feature_name, header, wk_dir, bam_file)
        out_file_list.append(data_file)