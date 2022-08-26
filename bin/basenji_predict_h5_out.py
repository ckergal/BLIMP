#!/usr/bin/env python

from __future__ import print_function

from optparse import OptionParser
from math import *

import json
import os

import h5py
import numpy as np
import pandas as pd
import pybedtools as pbd
import pyBigWig as pbw
import tensorflow as tf

if tf.__version__[0] == '1':
  tf.compat.v1.enable_eager_execution()

from corgi import bed
from corgi import dna_io
from corgi import seqnn_novar
from corgi import stream
from corgi import dataset

'''
basenji_predict_h5_out.py
'''

################################################################################
# main
################################################################################
def main():
  usage = 'usage: %prog [options] <params_file> <model_file> <bed_file>'
  parser = OptionParser(usage)
  parser.add_option('-f', dest='genome_fasta',
      default=None,
      help='Genome FASTA for sequences [Default: %default]')
  parser.add_option('-o', dest='out_dir',
      default='pred_out',
      help='Output directory [Default: %default]')
  parser.add_option('-p', dest='processes',
      default=None, type='int',
      help='Number of processes, passed by multi script')
  parser.add_option('--shifts', dest='shifts',
      default='0',
      help='Ensemble prediction shifts [Default: %default]')
  parser.add_option('-t', dest='targets_file',
      default=None, type='str',
      help='File specifying target indexes and labels in table format')
  parser.add_option('--expr', dest = 'targets_dir', default=None)
  parser.add_option('--length', dest = 'seq_length', type="int")
  parser.add_option('-s', dest = 'summarize', default = False,
    action='store_true',
    help='Summarize preds and expr into 1 coefficient')
  parser.add_option('-c', dest="crop",
    default = 0, type="int")
  parser.add_option('-i', dest="index",
    default = None, type="str")
  (options, args) = parser.parse_args()

  if len(args) == 3:
    params_file = args[0]
    model_file = args[1]
    bed_file = args[2]

  elif len(args) == 4:
    params_file = args[0]
    model_file = args[1]
    bed_file = args[2]
    ortho_bed = args[3]

  else:
    parser.error('Must provide parameter and model files and BED file')

  if not os.path.isdir(options.out_dir):
    os.mkdir(options.out_dir)

  options.shifts = [int(shift) for shift in options.shifts.split(',')]

  # update bed file
  if options.targets_dir :
    bed_pred_seq = pbd.BedTool(ortho_bed)
    bed_tar_seq = pbd.BedTool('%s/sequences.bed' % options.targets_dir) ## MODIF
    # bed_tar_seq = pbd.BedTool('%s/sequences_notrain.bed' % options.targets_dir) ## MODIF
    bed_pred_seq = pbd.BedTool(bed_pred_seq.intersect(bed_tar_seq, c = True))

    # print(bed_pred_seq, "   bed_pred_seq 1")

    pos_overlap = [elem[3] for elem in bed_pred_seq]
    # pos_overlap = [i for i,x in enumerate(pos_overlap) if x == '1']
    pos_overlap = [i for i,x in enumerate(pos_overlap)]
    # print(pos_overlap, "   pos_overlap \n")

    newbed = bed_file.split("/")[0:-1]
    newbed = "/".join(newbed)
    newbed = "%s/newbed%s.bed" % (newbed , options.index)

    newortho = ortho_bed.split("/")[0:-1]
    newortho = "/".join(newortho)
    newortho = "%s/newortho%s.bed" % (newortho, options.index)


    with open(bed_file, "r+") as f:
      old = f.readlines()
      new = [old[i] for i in pos_overlap]

    with open(newbed, "w") as f:
      for elem in new:
        f.write(elem)

    with open(ortho_bed, "r+") as f:
      old = f.readlines()
      new = [old[i] for i in pos_overlap]

    with open(newortho, "w") as f:
      for elem in new:
        f.write(elem)

    bed_file = newbed
    ortho_bed = newortho

    ortho_bed = pbd.BedTool(ortho_bed)
    bed_pred_seq = pbd.BedTool(ortho_bed.intersect(bed_tar_seq, wb = True))    


  #################################################################
  # read parameters and collect target information

  with open(params_file) as params_open:
    params = json.load(params_open)
  params_model = params['model']

  crop = options.crop

  if options.targets_file is None:
    target_slice = None
  else:
    targets_df = pd.read_table(options.targets_file, index_col=0)
    target_slice = targets_df.index

  #################################################################
  # setup model

  seqnn_model = seqnn_novar.SeqNN(params_model)
  seqnn_model.restore(model_file)
  seqnn_model.build_slice(target_slice, seq_length = options.seq_length)
  seqnn_model.build_ensemble(False, options.shifts)

  _, preds_length, preds_depth = seqnn_model.model.output.shape

  if type(preds_length) == tf.compat.v1.Dimension:
    preds_length = preds_length.value
    preds_depth = preds_depth.value


  preds_window = seqnn_model.model_strides[0]
  seq_crop = seqnn_model.target_crops[0]*preds_window


  #################################################################
  # construct model sequences. MODIF

  model_seqs_dna, model_seqs_coords = bed.make_bed_seqs(
    newortho, options.genome_fasta,
    params_model['seq_length'], stranded=False)

  num_seqs = len(model_seqs_dna)  

  # construct site coordinates
  site_seqs_coords = bed.read_bed_coords(bed_file, 128)  

  #################################################################
  # setup output

  assert(preds_length % 2 == 0)
  preds_mid = preds_length // 2
  
  if options.crop :
    seq_length = options.seq_length - (options.crop * 128 * 2)
  else :
    seq_length = options.seq_length

  site_preds_length = seq_length // preds_window
  site_preds_length = 8

  assert(site_preds_length % 2 == 0)
  site_preds_start = preds_mid - site_preds_length//2
  site_preds_end = site_preds_start + site_preds_length

  # initialize HDF5
  out_h5_file = '%s/predict%s.h5' % (options.out_dir, options.index)
  if os.path.isfile(out_h5_file):
    os.remove(out_h5_file)
  out_h5 = h5py.File(out_h5_file, 'w')


  # store site coordinates
  site_seqs_chr, site_seqs_start, site_seqs_end = zip(*site_seqs_coords)
  site_seqs_chr = np.array(site_seqs_chr)  
  site_seqs_start = np.array(site_seqs_start)
  site_seqs_end = np.array(site_seqs_end)
  out_h5.create_dataset('chrom', data=np.array(site_seqs_chr, dtype='S'))
  out_h5.create_dataset('start', data=site_seqs_start)
  out_h5.create_dataset('end', data=site_seqs_end)

  #################################################################
  # predict scores, write output

  dt = h5py.vlen_dtype(np.dtype('float16'))

  if options.summarize :
    out_h5.create_dataset(
      'preds',
      shape=(num_seqs, 1, preds_depth),
      dtype='float16')

  else : 
    out_h5.create_dataset(
      'preds',
      shape=(num_seqs, site_preds_length, preds_depth),
      dtype='float16')


  def seqs_gen():
    for seq_dna in model_seqs_dna :
      yield dna_io.dna_1hot(seq_dna)

  preds_stream = stream.PredStreamGen(seqnn_model, seqs_gen(), params['train']['batch_size'])


  # predict
  for si in range(num_seqs):
    preds_seq = preds_stream[si]
    
    site_preds_start = 508
    site_preds_end = 516

    preds_site = preds_seq[site_preds_start:site_preds_end,:]

    if options.summarize :
      out_h5['preds'][si] = sum(preds_site)

    else :
      out_h5['preds'][si] = preds_site


  #################################################################  MORE NEW THAN NEW
  
  site_seqs_coords = bed.read_bed_coords(newortho, 128) ## REVOIR
  site_seqs_chr, site_seqs_start, site_seqs_end = zip(*site_seqs_coords)

  if options.targets_dir :
    if options.summarize :
      out_h5.create_dataset(
        'expr',
        shape=(num_seqs, 1, preds_depth),
        dtype='float16')

    else :
      out_h5.create_dataset(
        'expr',
        shape=(num_seqs, site_preds_length, preds_depth),
        dtype='float16')

  toto = pd.read_table("%s/targets.txt" % options.targets_dir)
  for seq in range(num_seqs):
    for cov_file in range(len(toto['file'])) :
      cov_open = pbw.open(toto['file'][cov_file], 'r')    
      cov_values = cov_open.values(
        str(site_seqs_chr[seq]), 
        int(site_seqs_start[seq]) + site_preds_start * 128, 
        int(site_seqs_start[seq]) + site_preds_end * 128)
      expr = [sum(cov_values[i:i+128]) for i in range(0, len(cov_values), 128)]
      cov_open.close()

      if options.summarize: 
        out_h5['expr'][seq, :, cov_file] = sum(expr)

      else : 
        out_h5['expr'][seq, :, cov_file] = expr


  #################################################################  NEW
  # evluate

  # if options.targets_dir :
  #   if options.summarize :
  #     out_h5.create_dataset(
  #       'expr',
  #       shape=(num_seqs, 1, preds_depth),
  #       dtype='float16')

  #   else :
  #     out_h5.create_dataset(
  #       'expr',
  #       shape=(num_seqs, site_preds_length, preds_depth),
  #       dtype='float16')

  #   df_pred_seq = pd.DataFrame(
  #     {
  #     "chrPred" : [i[0] for i in bed_pred_seq],
  #     "startPred" : [i[1] for i in bed_pred_seq],
  #     "stopPred" : [i[2] for i in bed_pred_seq],

  #     "chrTar" : [i[3] for i in bed_pred_seq],
  #     "startTar" : [i[4] for i in bed_pred_seq],
  #     "stopTar" : [i[5] for i in bed_pred_seq]
  #     }
  #     )

  #   ### df_tar_seq : "whole genome" in h5 cov files
  #   df_tar_seq = pd.read_csv(
  #     '%s/sequences.bed' % options.targets_dir,
  #     sep="\t",
  #     usecols = [0,1,2],
  #     names = ["chr", "start", "stop"],
  #     header = None)

  #   target_slice = pd.read_csv(
  #     options.targets_file,
  #     sep = "\t")

  #   for cov_file in range(len(target_slice)) :
  #     tar_cov_h5 = h5py.File("%s/seqs_cov/%s.h5" % (options.targets_dir, cov_file), 'r') 
  #     tar_cov = pd.DataFrame(tar_cov_h5['targets'])
  #     tar_cov_h5.close()
  #     tar_cov = pd.concat([tar_cov, df_tar_seq], axis = 1)     

  #     count_seq = 0

  #     for seq in range(len(df_pred_seq)) :       

  #       if df_pred_seq.iloc[seq].startPred == df_pred_seq.iloc[seq - 1].stopPred :
  #         int_inf = (int(df_pred_seq.iloc[seq - 1].startPred) - int(df_pred_seq.iloc[seq - 1].startTar)) // 128
  #         int_sup = int_inf + ((int(df_pred_seq.iloc[seq].stopPred) - int(df_pred_seq.iloc[seq - 1].startPred)) // 128)

  #         int_inf += crop
  #         int_sup -= crop

  #         expr[-1] = list(tar_cov.loc[
  #           (tar_cov.start == int(df_pred_seq.iloc[seq].startTar)) &
  #           (tar_cov.chr == str(df_pred_seq.iloc[seq].chrTar))
  #           , list(range(int_inf,int_sup + 1))].values[0])          

  #         count_seq += 1

  #       else :

  #         int_inf = ((int(df_pred_seq.iloc[seq].startPred) - int(df_pred_seq.iloc[seq].startTar)) // 128)
  #         int_sup = (int_inf + ((int(df_pred_seq.iloc[seq].stopPred) - int(df_pred_seq.iloc[seq].startPred)) // 128))

  #         int_inf += crop
  #         int_sup -= crop

  #         real_seq = seq - count_seq

  #         list_pred = out_h5['preds'][real_seq, :, cov_file].tolist()

  #         expr = tar_cov.loc[
  #         (tar_cov.start == int(df_pred_seq.iloc[seq].startTar)) &
  #         (tar_cov.chr == str(df_pred_seq.iloc[seq].chrTar))
  #         , list(range(int_inf,int_sup))].values[0]          

  #         if options.summarize :             
  #           if len(list_pred) == len([sum(expr)]) and expr != [] :
  #             out_h5['expr'][seq, :, cov_file] = sum(expr)
  #           else : 
  #             out_h5['expr'][seq, :, cov_file] = np.nan
          
  #         else : 
  #           if len(list_pred) == len(expr) :
  #             out_h5['expr'][seq, :, cov_file] = expr
  #           else :
  #             out_h5['expr'][seq, :, cov_file] = np.nan


  #         # if len(list_pred) == len([sum(expr)]) and expr != [] :
  #         #   if options.summarize :
  #         #     out_h5['expr'][seq, :, cov_file] = sum(expr)
  #         #   else :
  #         #     out_h5['expr'][seq, :, cov_file] = expr

  #         # else :
  #         #   out_h5['expr'][seq, :, cov_file] = np.nan

  # out_h5.close()

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
  main()
