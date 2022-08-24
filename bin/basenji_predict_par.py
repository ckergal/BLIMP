#!/usr/bin/env python

from optparse import OptionParser

import os
import glob
import pybedtools as pbd
import h5py
import numpy as np
import pandas as pd

from scipy.stats import spearmanr, pearsonr, shapiro

from corgi import slurm

'''
basenji_predict_par.py
Predict variable sequences from a BED file.
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
    default=1, type='int',
    help='Number of processes, passed by multi script')
  parser.add_option('--shifts', dest='shifts',
    default='0',
    help='Ensemble prediction shifts [Default: %default]')
  parser.add_option('-t', dest='targets_file',
    default=None, type='str',
    help='File specifying target indexes and labels in table format')
  parser.add_option('--expr', dest='targets_dir', default=None)
  parser.add_option('--length', dest = 'seq_length', 
    type="int")
  parser.add_option('-s', dest = 'summarize', default = False,
    help='Summarize preds and expr into 1 coefficient', action='store_true')
  parser.add_option('-n', dest="acc_name",
    default="acc", type="str")
  parser.add_option('-c', dest="crop",
    default=0, type="int")
  parser.add_option('-l', dest="lines",
    default=100, type="int")
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

  ortho_bed = pbd.BedTool(ortho_bed)
  bed_file = pbd.BedTool(bed_file)

  steps = [i for i in range(0, len(ortho_bed), options.lines)]
  steps.append(len(ortho_bed))  

  ### modif pour test
  # bed_file = ortho_bed

  sequences=[]
  sequences_ortho=[]
  if len(bed_file) >= options.lines : 
    for i in range(len(steps)-1):
      sequences.append([str(seq) for seq in bed_file[steps[i]:steps[i+1]]])
      sequences_ortho.append([str(seq) for seq in ortho_bed[steps[i]:steps[i+1]]])
  
  for i in range(len(sequences)):
    with open("%s/%s.bed" % (options.out_dir, i), "w") as f:
      for seq in sequences[i]:
        f.write(seq)

  for i in range(len(sequences_ortho)):
    with open("%s/%s_ortho.bed" % (options.out_dir, i), "w") as f:
      for seq in sequences_ortho[i]:
        f.write(seq)



  ortho_files = glob.glob("%s/*ortho.bed" % options.out_dir)
  seq_files = glob.glob("%s/[0-9A-Fa-f].bed" % options.out_dir)

  # print(seq_files)
  # print(ortho_files)
  # exit()
  read_jobs = []  

  for index_file in range(len(seq_files)) :
    cmd = 'basenji_predict_h5_out.py'
    cmd += ' -o %s' % options.out_dir
    cmd += ' -t %s' % options.targets_file
    cmd += ' -f %s' % options.genome_fasta
    cmd += ' --expr %s' % options.targets_dir
    cmd += ' -c %s' % options.crop
    cmd += ' -i %d' % index_file
    cmd += ' --length %d' % options.seq_length
    if options.summarize :
      cmd += ' -s '
    cmd += ' %s' % params_file
    cmd += ' %s' % model_file
    cmd += ' %s/%s.bed' % (options.out_dir, index_file)
    cmd += ' %s/%s_ortho.bed' % (options.out_dir, index_file)

    j = slurm.Job(cmd,
      name='seq%s' % index_file,
      out_file='%s/%s.out' % (options.out_dir, index_file),
      #err_file='%s/%s.err' % (options.out_dir, index_file),
      queue='standard',
      # mem=15000,
      time='12:0:0')
    read_jobs.append(j)

  slurm.multi_run(read_jobs, options.processes, verbose = True,
    launch_sleep=1, update_sleep=5)

  # exit()

  # ################################################################
  # # read preds and expr
  seqs_cov_files = []
  ti = 0
  seqs_cov_file = '%s/predict%d.h5' % (options.out_dir, ti)

  while os.path.isfile(seqs_cov_file):
    seqs_cov_files.append(seqs_cov_file)
    ti += 1    
    seqs_cov_file = '%s/predict%d.h5' % (options.out_dir, ti)

  os.system("cat %s/newortho*.bed > %s/ortho_species.bed" %(options.out_dir, options.out_dir))
  os.system("cat %s/newbed*.bed > %s/model_species.bed" %(options.out_dir, options.out_dir))
  os.system("rm %s/new*.bed" %options.out_dir)

  ortho_bed = pbd.BedTool("%s/ortho_species.bed" % options.out_dir)


  ################################################################
  # correlation by sample  
  corr = []
  samples_file = pd.read_table(options.targets_file, index_col = 0)

  for sample_index in range(samples_file.shape[0]) :
    pred = []
    expr = []
    for file in seqs_cov_files :
      res_open = h5py.File(file, 'r')
      pred.append(res_open['preds'][:,:,sample_index].flatten())      
      expr.append(res_open['expr'][:,:,sample_index].flatten())      

      res_open.close()

    expr = np.concatenate(expr, axis=None)
    pred = np.concatenate(pred, axis=None)

    pred = pred[np.logical_not(np.isnan(expr))]
    expr = expr[np.logical_not(np.isnan(expr))]

    

    if shapiro(pred)[1] > .05 and shapiro(expr)[1] > .05 :
      corr.append(pearsonr(pred, expr)[0])

    else:
      corr.append(spearmanr(expr, pred)[0])


  corr_sample = pd.DataFrame({
    'index':samples_file.index,
    'corr':corr
    })

  corr_sample.to_csv("%s/corr_sample.txt" % options.out_dir,
    sep = "\t",
    index = False,
    float_format="%.5f")
  
  
  ################################################################
  # correlation by sequence
  corr = []
  for file in seqs_cov_files:
    res_open = h5py.File(file, 'r')

    for seq_index in range(res_open['preds'].shape[0]):

      pred = np.array(res_open['preds'][seq_index,:,:].flatten(), dtype='f8')      
      expr = np.array(res_open['expr'][seq_index,:,:].flatten(), dtype='f8')      


      if shapiro(pred)[1] > .05 and shapiro(expr)[1] > .05:
        corr.append(pearsonr(pred, expr)[0])

      else :
        corr.append(spearmanr(pred, expr)[0])

    res_open.close()

  corr_seq = pd.DataFrame({
    'chr': [seq[0] for seq in ortho_bed],
    'start': [seq[1] for seq in ortho_bed],
    'stop': [seq[2] for seq in ortho_bed],
    'corr': corr
    })

  corr_seq.to_csv("%s/corr_seq.txt" % options.out_dir,
    sep="\t",
    index=False,
    float_format='%.5f')


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
  main()