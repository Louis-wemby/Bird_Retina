import matplotlib.pyplot as plt
import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = pow(2,40).__str__()
import cv2
from scipy import ndimage as ndi
from skimage import (
    color, feature, filters, measure, segmentation, io
)
import numpy as np
import pandas as pd
import os


class Segmentation:

    def __init__(self):
        self.img_path = None
        self.mRNA_path = None

        self.raw_img = None
        self.img = None
        self.mask = None
        self.label = None

    def load(self, img_path, mRNA_path, signal_pbar=None):
        self.mRNA_path = mRNA_path
        self.img_path = img_path
        self.raw_img = cv2.imread(self.img_path, cv2.IMREAD_GRAYSCALE)

    def pre_process(self,
                    threshold='auto',
                    verbose=True,
                    signal_pbar=None
                    ):
        if threshold == 'auto':
            threshold, _ = cv2.threshold(self.raw_img.copy(), 0, 255, cv2.THRESH_OTSU)
        if threshold > 0:
            _, self.img = cv2.threshold(self.raw_img.copy(), threshold, 255, cv2.THRESH_TOZERO)
        else:
            self.img = self.raw_img.copy()
        if verbose:
            print(f'Used Threshold: {threshold}')
            plt.figure(figsize=(16, 16))
            plt.imshow(self.img, 'gray')

    def watershed(self,
                  block_size=41,
                  offset=0.003,
                  min_distance=15,
                  expand_distance=0,
                  verbose=True,
                  signal_pbar=None
                  ):
        img = self.img.copy()
        threshold = filters.threshold_local(img, block_size=block_size, offset=offset)
        if verbose:
            plt.figure(figsize=(16, 16))
            plt.imshow(threshold)
            plt.title('Local threshold')
        distance = ndi.distance_transform_edt(img > threshold)
        if verbose:
            plt.figure(figsize=(16, 16))
            plt.imshow(distance)
            plt.title('Distance map')

        local_max_coords = feature.peak_local_max(distance, min_distance=min_distance)
        local_max_mask = np.zeros(distance.shape, dtype=bool)
        local_max_mask[tuple(local_max_coords.T)] = True
        markers = measure.label(local_max_mask)
        if verbose:
            plt.figure(figsize=(16, 16))
            plt.imshow(markers)
            plt.title('Markers')

        self.mask = segmentation.watershed(-distance, markers, mask=img)
        if expand_distance > 0:
            self.mask = segmentation.expand_labels(self.mask, distance=expand_distance)
        label = color.label2rgb(self.mask, bg_label=0)
        self.label = (label * 255).astype(np.uint8)
        if verbose:
            print("numbers of cells:", self.mask.max())
            plt.figure(figsize=(16, 16))
            plt.imshow(self.label)
            plt.title('Label image')

    def cellpose(self,
                 gpu=False,
                 pretrained_model='cyto',
                 diameter=None,
                 flow_threshold=0,
                 mask_threshold=-2,
                 verbose=True,
                 expand_distance=0,
                 min_size=200,  #
                 signal_pbar=None,
                 **kwargs
                 ):
        from cellpose import models, plot, utils
        # define available model names, right now we have three broad categories
        model_names = ['cyto', 'nuclei', 'bact', 'cyto2', 'bact_omni', 'cyto2_omni']
        builtin_model = np.any([pretrained_model == s for s in model_names])

        img = self.img.copy()
        if builtin_model:
            model = models.Cellpose(gpu=gpu, model_type=pretrained_model)
        else:
            cpmodel_path = pretrained_model
            model = models.CellposeModel(gpu=gpu, pretrained_model=cpmodel_path)
        out = model.eval([img], diameter=diameter, channels=[0, 0],
                         flow_threshold=flow_threshold, mask_threshold=mask_threshold,
                         do_3D=False, **kwargs)
        masks, flows = out[:2]

        if len(out) > 3:
            diams = out[-1]

        self.mask = masks[0]
        if expand_distance > 0:
            self.mask = segmentation.expand_labels(self.mask, distance=expand_distance)
        if min_size > 0:
            self.mask = utils.fill_holes_and_remove_small_masks(self.mask, min_size=min_size)
        label = color.label2rgb(self.mask, bg_label=0)
        self.label = (label * 255).astype(np.uint8)

        outlines = utils.masks_to_outlines(self.mask)
        outX, outY = np.nonzero(outlines)
        img0 = np.zeros((self.img.shape[0], self.img.shape[1], 3), dtype=np.uint8)
        img0[:, :, 0] = self.img
        img0[:, :, 1] = self.img
        img0[:, :, 2] = self.img

        img0[outX, outY, :] = np.array([255, 0, 0])  # pure red
        self.img_lines = img0

        if verbose:
            print("numbers of cells:", self.mask.max())
            plt.figure(figsize=(16, 16))
            plt.imshow(self.label)
            plt.title('Labels')

            plt.figure(figsize=(16, 16))
            plt.imshow(self.img_lines)
            plt.title('Outlines')

    def save_scGEM(self,
                   save_path,
                   name,
                   verbose=True,
                   signal_pbar=None,
                   minus_min=True, 
                   ):
        data = pd.read_csv(self.mRNA_path, sep='\t', comment="#")

        seg_cell_coor = []
        min_x = data['x'].min() if minus_min else 0
        min_y = data['y'].min() if minus_min else 0
        for i in range(self.mask.shape[0]):
            for j in range(self.mask.shape[1]):
                c = self.mask[i, j]
                if c:
                    seg_cell_coor.append([i + min_x, j + min_y, c])
        if signal_pbar:
            signal_pbar.emit(70)
        seg_cell_coor = pd.DataFrame(seg_cell_coor, columns=['x', 'y', 'cell'])
        cell_data = pd.merge(data, seg_cell_coor, how='left', on=['x', 'y'])
        cell_data = cell_data.dropna()
        cell_data['cell'] = cell_data['cell'].astype(int)
        # name = os.path.basename(self.mRNA_path)
        # name = os.path.splitext(name)[0]
        mask_fn = os.path.join(save_path, f'{name}_mask.npy')
        np.save(mask_fn, self.mask)
        gem_fn = os.path.join(save_path, f'{name}_scgem.csv.gz')
        cell_data.to_csv(gem_fn, index=False, sep='\t', compression="gzip")
        # coor_fn = os.path.join(save_path, f'{name}.ssDNA_coor.csv')
        # seg_cell_coor.to_csv(os.path.join(save_path, f'{args.i}.ssDNA_coor.csv'), index=False)
        if verbose:
            print(f'segmented mask save path: {mask_fn}')
            print(f'single-cell GEM save path: {gem_fn}')

    def __repr__(self):
        t = f"ssDNA Image Segmentation Object\n" \
            f"Raw   Image Path: {self.img_path}\n" \
            f"GEM   Data  Path: {self.mRNA_path}"
        return t
