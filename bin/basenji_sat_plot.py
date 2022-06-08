#!/usr/bin/env python
# Copyright 2017 Calico LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# =========================================================================
from __future__ import print_function

from optparse import OptionParser
import os
import pdb
import random
import sys

import h5py
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from PIL import Image
import seaborn as sns

'''
basenji_sat_plot.py
Generate plots from scores HDF5 file output by saturation mutagenesis analysis
via basenji_sat_bed.py
'''

################################################################################
# main
################################################################################
def main():
  usage = 'usage: %prog [options] <scores_file>'
  parser = OptionParser(usage)
  parser.add_option('-f', dest='figure_width',
      default=20, type='float',
      help='Figure width [Default: %default]')
  parser.add_option('-l', dest='plot_len',
      default=200, type='int',
      help='Length of centered sequence to mutate [Default: %default]')
  parser.add_option('-m', dest='min_limit',
      default=0.05, type='float',
      help='Minimum heatmap limit [Default: %default]')
  parser.add_option('-o', dest='out_dir',
      default='sat_plot', help='Output directory [Default: %default]')
  parser.add_option('-r', dest='rng_seed',
      default=1, type='float',
      help='Random number generator seed [Default: %default]')
  parser.add_option('-s', dest='sample',
      default=None, type='int',
      help='Sample N sequences from the set [Default:%default]')
  parser.add_option('--stat', dest='sad_stat',
      default='sum',
      help='SAD stat to display [Default: %default]')
  parser.add_option('-t', dest='targets_file',
      default=None, type='str',
      help='File specifying target indexes and labels in table format')
  parser.add_option('--ts', dest='target_sample',
    default=None, help='Coma-separated list of target indexes to make plots')
  parser.add_option('--tm', dest='target_muts',
    default = None, help='Coma-separated known mutations targeted')  
  (options, args) = parser.parse_args()

  scores_h5_file = args[0]

  if not os.path.isdir(options.out_dir):
    os.mkdir(options.out_dir)

  save_ext = 'png'

  np.random.seed(options.rng_seed)

  # open scores
  scores_h5 = h5py.File(scores_h5_file, mode = 'r')


  # extract shapes
  num_seqs = scores_h5['seqs'].shape[0]
  mut_len = scores_h5[options.sad_stat].shape[1]

  if options.plot_len > mut_len:
    print('Decreasing plot_len=%d to maximum %d' % (options.plot_len, mut_len), file=sys.stderr)
    options.plot_len = mut_len

  # determine targets
  if options.targets_file is not None:
    targets_df = pd.read_table(options.targets_file, index_col=0)
    num_targets = targets_df.shape[0]
  else:
    num_targets = scores_h5[options.sad_stat].shape[-1]


  # determine plot region
  mut_mid = mut_len // 2
  plot_start = mut_mid - (options.plot_len//2)
  plot_end = plot_start + options.plot_len

  # plot attributes
  sns.set(style='white', font_scale=1)
  spp = subplot_params(options.plot_len)

  # determine sequences
  seq_indexes = np.arange(num_seqs)

  if options.sample and options.sample < num_seqs:
    seq_indexes = np.random.choice(seq_indexes, size=options.sample, replace=False)

  if options.target_sample is not None:
    target_sample = [int(ti) for ti in options.target_sample.split(',')]

  if options.target_muts is not None : 
    target_muts = [int(ti) for ti in options.target_muts.split(',')]
    
    for si in seq_indexes:
      index_muts = []
      for mut in target_muts :
        index_muts.append(mut - scores_h5['start'][si] - 1)
      
      # read sequence
      seq_1hot = scores_h5['seqs'][si,plot_start:plot_end]

      # read scores
      scores = scores_h5[options.sad_stat][si,plot_start:plot_end,:,:]

      # reference scores
      ref_scores = scores[seq_1hot]

      for ti in range(num_targets):
        
        scores_ti = scores[:,:,ti]

        # compute scores relative to reference
        delta_ti = scores_ti - ref_scores[:,[ti]]
  
        # compute loss and gain
        delta_loss = delta_ti.min(axis=1)
        delta_gain = delta_ti.max(axis=1)

        # setup plot
        fig = plt.figure(figsize=(options.figure_width, 6))
        fig.suptitle(targets_df.iloc[ti]['description'], fontsize=16)
        
        grid_rows = 2
        row_i = 0        
        ax_sad = plt.subplot2grid(
          (grid_rows, spp['heat_cols']), (row_i, spp['sad_start']),
          colspan=spp['sad_span'])
        row_i += 1
        ax_heat = plt.subplot2grid(
          (grid_rows, spp['heat_cols']), (row_i, 0), colspan=spp['heat_cols'])
     
        # plot SAD
        plot_sad(ax_sad, delta_loss, delta_gain, target_muts=target_muts, index_muts = index_muts)

        # plot heat map
        plot_heat(ax_heat, delta_ti.T, options.min_limit, target_muts=index_muts)

        plt.tight_layout()        
        plt.savefig('%s/seq%d_%s.%s' % (options.out_dir, si, targets_df.index[ti], save_ext), dpi=100)
        plt.close()




def plot_heat(ax, sat_delta_ti, min_limit, target_muts = None):
  """ Plot satmut deltas.
    Args:
        ax (Axis): matplotlib axis to plot to.
        sat_delta_ti (4 x L_sm array): Single target delta matrix for saturated mutagenesis region,
        min_limit (float): Minimum heatmap limit.
    """    
  for mut in target_muts:
    plt.axvline(x=mut, c="black", ls=":")

  vlim = max(min_limit, np.nanmax(np.abs(sat_delta_ti)))
  sns.heatmap(
      sat_delta_ti,
      linewidths=0,
      cmap='RdBu_r',
      vmin=-vlim,
      vmax=vlim,
      xticklabels=False,
      ax=ax)
  ax.yaxis.set_ticklabels('ACGT', rotation='horizontal')  # , size=10)



def plot_sad(ax, sat_loss_ti, sat_gain_ti, target_muts=None, index_muts = None):
  """ Plot loss and gain SAD scores.
    Args:
        ax (Axis): matplotlib axis to plot to.
        sat_loss_ti (L_sm array): Minimum mutation delta across satmut length.
        sat_gain_ti (L_sm array): Maximum mutation delta across satmut length.
    """

  rdbu = sns.color_palette('RdBu_r', 10)

  ax.plot(-sat_loss_ti, c=rdbu[0], label='loss', linewidth=1)
  ax.plot(sat_gain_ti, c=rdbu[-1], label='gain', linewidth=1)
  ax.set_xlim(0, len(sat_loss_ti))
  ax.legend()  

  for mut in index_muts:
    ax.axvline(x=mut, c="black", ls=":")

  for i in range(len(index_muts)):
    ax.text(index_muts[i], 1.1*max(sat_gain_ti), str(target_muts[i]),fontweight="bold", fontsize=8)    

  ax.xaxis.set_ticks([])
  for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(0.5)


def subplot_params(seq_len):
  """ Specify subplot layout parameters for various sequence lengths. """
  if seq_len < 500:
    spp = {
        'heat_cols': 400,
        'sad_start': 1,
        'sad_span': 321,
        'logo_start': 0,
        'logo_span': 323
    }
  else:
    spp = {
        'heat_cols': 400,
        'sad_start': 1,
        'sad_span': 320,
        'logo_start': 0,
        'logo_span': 322
    }

  return spp


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
  main()
