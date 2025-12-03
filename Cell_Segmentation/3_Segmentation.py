import sys
sys.path.append('./RYTools')
from Tools.Segmentation import Segmentation
import matplotlib.pyplot as plt
from Tools.Spatial import *
import numpy as np
import pandas as pd
import scanpy as sc
import os
import warnings
warnings.filterwarnings('ignore')

inputfile = sys.argv[1]
outputfile = sys.argv[2]
idname = sys.argv[3]
bs = sys.argv[4]
ot = sys.argv[5]
md = sys.argv[6]
ed = sys.argv[7]

bs = int(bs)  # block size
ot = float(ot)  # offset
md = int(md)  # minimum distance
ed = int(ed)  # expanding distance

if not os.path.exists(outputfile):
    os.mkdir(outputfile)

# initializing segmentation object
sobj = Segmentation()
sobj.load(
    img_path = f'{inputfile}/{idname}/{idname}_matched_ssDNA.png', 
    mRNA_path = f'{inputfile}/{idname}/{idname}.gem.gz',
)

plt.figure(figsize=(16,16))
plt.imshow(sobj.raw_img, 'gray')

# pre processing   sobj.raw_img --> sobj.img
sobj.pre_process(threshold='auto')

# watershed algorithm
sobj.watershed(
    block_size=bs,
    offset=ot,
    min_distance=md,
    expand_distance=ed,
    verbose=False
)

result_path = f'{outputfile}/{idname}_label_ssDNA.png'
cv2.imwrite(result_path, sobj.label)

# save as scgem.csv format
save_path = os.path.join(outputfile)
sobj.save_scGEM(save_path = save_path, name = idname)
