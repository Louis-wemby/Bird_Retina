import numpy as np
import pandas as pd
import cv2
import sys
import os

idname = sys.argv[1]
inputdir = sys.argv[2]
outputdir = sys.argv[3]

if not os.path.exists(idname):
    os.mkdir(idname)

data = pd.read_csv(os.path.join(inputdir, f'{idname}.tissue.gem'), sep='\t', comment='#')
data['x'] -= data['x'].min()
data['y'] -= data['y'].min()
data = data.groupby(['x', 'y']).sum().reset_index()
w = data['x'].max() + 1
h = data['y'].max() + 1
mtx = np.zeros((w, h))
mtx[data['x'], data['y']] = data['MIDCount']
mtx *= 25
mtx = np.clip(mtx, 0, 255)
cv2.imwrite(os.path.join(outputdir, f'{idname}.mRNA.png'), mtx)

