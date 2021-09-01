
# try saving the sparse matrix as mtx
import scipy.sparse as sparse
import scipy.io as sio
import scanpy as sc

import sys, getopt

def main(argv):
   inputfile = ''
   outputfile = ''
   layer_id = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:c:",["ifile=","ofile=","counts="])
   except getopt.GetoptError:
      print('test.py -i <inputfile> -o <outputfolder>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('test.py -i <inputfile> -o <outputfolder>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
      elif opt in ("-c", "--counts"):
         layer_id = arg

   print('Input file is "' + inputfile)
   print('Output folder is "' + outputfile)

   adata = sc.read_h5ad(inputfile)
   adata.write_csvs( outputfile +  "/", skip_data=True)
   # Select which matrix to export. Raw counts should be default behaviour
   # if layer_id == log then we use the log normalized data. 
   if(layer_id!='log'):
       adata.X = adata.layers[layer_id].copy()

   sio.mmwrite(outputfile +"/sparse_matrix.mtx",adata.X)


if __name__ == "__main__":
   main(sys.argv[1:])
