# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 18:22:47 2021
RY Spatial Analysis 
@author: wangshuai3
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from anndata import AnnData
from scipy.spatial import cKDTree
from typing import Optional, List, Tuple, Union
import base64
import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = pow(2,40).__str__()
os.environ["OPENCV_SHRT_MAX"] = str(pow(2,40))
import cv2


def read_loom(fname: str, npy: str, id_: str) -> sc.AnnData:
    """
    

    Parameters
    ----------
    fname : str
        loom filename
    npy : str
        npy filename
    id_ : str
        id

    Returns
    -------
    obj : TYPE
        AnnData

    """
    obj = sc.read(fname)
    obj.uns['id_'] = id_
    obj.uns['seg_cell'] = np.load(npy)
    obj.obs['cell_id'] = obj.obs_names.to_series()
    obj.obs['cell_id'] = obj.obs['cell_id'].apply(lambda x: int(x.split('_')[1]))
    obj.obs['cell_id'] = pd.Categorical(obj.obs['cell_id'])
    obj = get_cell_center(obj)
    return obj


def black_line(arr, arr2, col=(0, 0, 0, 1), sep=2):
    # arr += 1
    arr_pad = np.pad(arr, sep, constant_values=-1)
    arr_up = np.pad(arr, [(0, sep * 2), (sep, sep * 1)], constant_values=-1)
    arr_left = np.pad(arr, [(sep, sep), (0, sep * 2)], constant_values=-1)
    arr_down = np.pad(arr, [(sep * 2, 0), (sep, sep * 1)], constant_values=-1)
    arr_right = np.pad(arr, [(sep, sep), (sep * 2, 0)], constant_values=-1)
    for i in [arr_up, arr_left, arr_down, arr_right]:
        arr2[(i - arr_pad)[sep:-sep, sep:-sep] != 0] = col
    # arr -= 1
    return arr2


def Hex_to_RGB(hex):
    r = int(hex[1:3], 16)
    g = int(hex[3:5], 16)
    b = int(hex[5:7], 16)
    return r, g, b


colorlist = ["#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941",
             "#006FA6", "#A30059", "#FFE4E1", "#0000A6", "#63FFAC",
             "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007",
             "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#FF2F80",
             "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9",
             "#B903AA", "#DDEFFF", "#7B4F4B", "#A1C299", "#0AA6D8",
             "#00A087FF", "#4DBBD5FF", "#E64B35FF", "#3C5488FF", "#F38400",
             "#A1CAF1", "#C2B280", "#848482", "#E68FAC", "#0067A5",
             "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300",
             "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26",
             "#191970", "#000080",
             "#6495ED", "#1E90FF", "#00BFFF", "#00FFFF", "#FF1493",
             "#FF00FF", "#A020F0", "#63B8FF", "#008B8B", "#54FF9F",
             "#00FF00", "#76EE00", "#FFF68F"]

colorlist = [Hex_to_RGB(i) for i in colorlist]


def cut_mtx(mtx, sep=0, constant_values=-1, threshold=0, axis=2):
    x1, x2 = np.where(mtx.sum(axis=axis) > threshold)
    mtx = mtx[x1.min():x1.max() + 1, x2.min():x2.max() + 1].copy()
    mtx = np.pad(mtx, sep, mode='constant', constant_values=constant_values)
    return mtx


def save_svg(fname, arr, ah, aw):
    plt.imsave(fname + '.png', arr)
    startSvgTag = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"
"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg version="1.1"
xmlns="http://www.w3.org/2000/svg"
xmlns:xlink="http://www.w3.org/1999/xlink"
"""
    figTag = f'width="{ah}px" height="{aw}px" viewBox="0 0 {ah} {aw}">'
    endSvgTag = """</svg>"""
    with open(fname + '.png', 'rb') as f:
        data = f.read()
    base64String = f'<image xlink:href="data:image/png;base64,{base64.b64encode(data).decode()}" width="{ah}" height="{aw}" x="0" y="0" />'
    with open(fname, 'w') as f:
        f.write(startSvgTag + figTag + base64String + endSvgTag)
    os.remove(fname + '.png')


def featureplot_cell_discrete(obj: sc.AnnData,
                              feature: str,
                              fname: Union[None, str] = None,
                              order: List = None,
                              colors: List = None,
                              bg_color=(0, 0, 0),
                              line_color=(0, 0, 0),
                              blank_color=(192, 192, 192),
                              show: bool = True,
                              title: Union[None, str] = None,
                              scale: bool = False,
                              sep: int = 2,
                              dpi: int = 600,
                              legend_size: int = 12
                              ):
    """
    

    Parameters
    ----------
    obj : sc.AnnData
        DESCRIPTION.
    feature : str
        selected feature ,colname or gene name of obs（obj.var_names）
    fname : TYPE, optional
        save filename (None: not save)
    order : TYPE, optional
        to draw in specific type & order
    colors : TYPE, optional
        select color for specific type & order to draw
    bg_color : TYPE, optional
        DESCRIPTION. The default is (0,0,0).
    line_color : TYPE, optional
        DESCRIPTION. The default is (0,0,0).
    blank_color : TYPE, optional
        DESCRIPTION. The default is (192,192,192).
    show : TYPE, optional
        DESCRIPTION. The default is True : bool.
    title : TYPE, optional
        DESCRIPTION. The default is None : Union[None, str].
    scale : TYPE, optional
        DESCRIPTION. The default is False : bool.
    sep : TYPE, optional
        DESCRIPTION. The default is 2 : int.
    dpi : TYPE, optional
        DESCRIPTION. The default is 600 : int.
    legend_size : TYPE, optional
        all cell types set to 6

    Returns
    -------
    None.

    """
    if colors is None:
        colors = colorlist
    elif isinstance(colors[0], str):
        colors = [Hex_to_RGB(i) for i in colors]
    data = obj.obs[feature].astype(str)
    clusters = data.unique().tolist()
    if order:
        if len(order) < len(clusters):
            clusters = order
            colors = colors[:len(clusters)]
            data[~data.isin(order)] = ''
            clusters.append('')
            colors.append(blank_color)
        else:
            clusters = order
            colors = colors[:len(clusters)]
            colors = colors[:len(clusters)]
    else:
        colors = colors[:len(clusters)]
    col_df = pd.DataFrame(dict(clusters=clusters,
                               r=[i[0] for i in colors],
                               g=[i[1] for i in colors],
                               b=[i[2] for i in colors]))
    data = pd.merge(data, col_df, left_on=feature, right_on='clusters', how='left')
    data['cell'] = obj.obs['cell_id'].values
    seg_df = pd.DataFrame(obj.uns['seg_cell'].flatten(), columns=['cell'])
    seg_df = pd.merge(seg_df, data, how='left', on='cell')
    seg_df.loc[seg_df['r'].isnull(), ['r', 'g', 'b']] = bg_color
    d1, d2 = obj.uns['seg_cell'].shape
    seg_display = seg_df[['r', 'g', 'b']].to_numpy().reshape(d1, d2, 3).astype(np.uint8)
    seg_display = black_line(obj.uns['seg_cell'], seg_display, sep=sep, col=line_color)

    fig = plt.figure(figsize=(16, 17), dpi=dpi)
    gs = gridspec.GridSpec(16, 17)
    if len(colors) < 16:
        ax1 = fig.add_subplot(gs[0, :len(colors)])
    else:
        ax1 = fig.add_subplot(gs[0, :])
    ax1.axis('off')
    ax1.imshow(np.array(colors).astype(np.uint8).reshape(1, len(colors), 3), aspect='auto')
    for ind, cluster in enumerate(clusters):
        ax1.text(ind - 0.5, -0.6, cluster.replace(' ', '\n'), fontdict={'size': str(legend_size)}, ha='left')
    ax2 = fig.add_subplot(gs[1:, :])
    ax2.imshow(seg_display)
    if not scale:
        ax2.axis('off')
    if title is None:
        title = feature
    ax2.set_title(title)
    if fname is not None:
        plt.savefig(fname)
        print('save as', fname)
    if not show:
        plt.close()


def featureplot_cell_continuous(obj,
                                feature,
                                fname=None,
                                cmap='viridis',
                                bg_color=(0.0, 0.0, 0.0, 1.0),
                                line_color=(0.0, 0.0, 0.0, 1.0),
                                show=True,
                                max_=None,
                                min_=None,
                                title=None,
                                scale=False,
                                sep=2,
                                dpi=600
                                ):
    """
    

    Parameters
    ----------
    obj : TYPE
        DESCRIPTION.
    feature : TYPE
        DESCRIPTION.
    fname : TYPE, optional
        DESCRIPTION. The default is None.
    cmap : TYPE, optional
        camp name，cmap in plt. The default is 'viridis'： str.
    bg_color : TYPE, optional
        added Channel A(0,0,0,0), transparent, float between 0-1. The default is (0,0,0,1).
    line_color : TYPE, optional
        DESCRIPTION. The default is (0,0,0,1).
    show : TYPE, optional
        DESCRIPTION. The default is True.
    max_ : TYPE, optional
        DESCRIPTION. The default is None.
    min_ : TYPE, optional
        DESCRIPTION. The default is None.
    title : TYPE, optional
        DESCRIPTION. The default is None.
    scale : TYPE, optional
        DESCRIPTION. The default is False.
    sep : TYPE, optional
        DESCRIPTION. The default is 2.
    dpi : TYPE, optional
        DESCRIPTION. The default is 600.

    Returns
    -------
    None.

    """
    cmap = plt.cm.get_cmap(cmap)
    if feature in obj.obs_keys():
        data = obj.obs.loc[:, [feature]]
    else:
        arr = obj[:, feature].X.toarray()
        arr = arr.reshape(arr.shape[0])
        data = pd.DataFrame(arr, index=obj.obs_names, columns=[feature])
    data['cell'] = obj.obs['cell_id'].values
    if max_ is not None:
        data.loc[data[feature] > max_, feature] = max_
    else:
        max_ = data[feature].max()
    if min_ is not None:
        data.loc[data[feature] < min_, feature] = min_
    else:
        min_ = data[feature].min()
    data[feature] = (data[feature] - min_) / max_
    data[feature] = data[feature].fillna(0)
    data[['r', 'g', 'b', 'a']] = cmap(data[feature])
    seg_df = pd.DataFrame(obj.uns['seg_cell'].flatten(), columns=['cell'])
    seg_df = pd.merge(seg_df, data, how='left', on='cell')
    seg_df.loc[seg_df['r'].isnull(), ['r', 'g', 'b', 'a']] = [bg_color]
    d1, d2 = obj.uns['seg_cell'].shape
    seg_display = seg_df[['r', 'g', 'b', 'a']].to_numpy().reshape(d1, d2, 4)
    seg_display = black_line(obj.uns['seg_cell'], seg_display, sep=sep, col=line_color)

    fig = plt.figure(figsize=(16, 17), dpi=dpi)
    gs = gridspec.GridSpec(16, 17)
    ax1 = fig.add_subplot(gs[0, :5])
    ax1.axis('off')
    ax1.imshow(cmap(np.linspace(0, 1, 100)).reshape(1, 100, 4), aspect='auto')
    ax1.text(-0.5, -0.6, f'{min_:.3f}', fontdict={'size': '12'}, ha='left')
    ax1.text(100, -0.6, f'{max_:.3f}', fontdict={'size': '12'}, ha='right')
    ax2 = fig.add_subplot(gs[1:, :])
    ax2.imshow(seg_display)
    if not scale:
        ax2.axis('off')
    if title is None:
        title = feature
    ax2.set_title(title)
    if fname is not None:
        plt.savefig(fname)
        print('save as', fname)
    if not show:
        plt.close()
    return max_, min_


def rotate_bound_white_bg(image, angle, bg=(0, 0, 0)):
    (h, w) = image.shape[:2]
    (cX, cY) = (w // 2, h // 2)
    M = cv2.getRotationMatrix2D((cX, cY), -angle, 1.0)
    cos = np.abs(M[0, 0])
    sin = np.abs(M[0, 1])
    nW = int((h * sin) + (w * cos))
    nH = int((h * cos) + (w * sin))
    M[0, 2] += (nW / 2) - cX
    M[1, 2] += (nH / 2) - cY
    return cv2.warpAffine(image, M, (nW, nH), borderValue=bg)


def featureplot_slices_discrete(obj,
                                feature,
                                fname=None,
                                order=None,
                                colors=None,
                                bg_color=(0, 0, 0),
                                line_color=(0, 0, 0),
                                blank_color=(192, 192, 192),
                                show=True,
                                scale=False,
                                title=None,
                                dpi=600,
                                sep=2,
                                fig_sep=0,
                                nrow=1,
                                ncol=1,
                                slices=None,
                                angle_dict={},
                                compress_factor=None,
                                raw=True,
                                legend_size=12
                                ):
    """
    angle_dict:
        obj.uns['angle_dict']
    raw:
        select when editing photo & visualizing web , no legend
    fame:
        web visualization ends with .svg
    """
    if slices is None:
        obj = obj.copy()
        slices = obj.obs['Batch'].unique().tolist()
    else:
        obj = obj[obj.obs['Batch'].isin(slices), :].copy()
    if colors is None:
        colors = colorlist
    elif isinstance(colors[0], str):
        colors = [Hex_to_RGB(i) for i in colors]
    data = obj.obs[[feature, 'Batch']]
    data[feature] = data[feature].astype(str)
    clusters = data[feature].unique().tolist()
    if order:
        if len(order) < len(clusters):
            clusters = order
            colors = colors[:len(clusters)]
            data.loc[~data[feature].isin(order), feature] = ''
            clusters.append('')
            colors.append(blank_color)
        else:
            clusters = order
            colors = colors[:len(clusters)]
            colors = colors[:len(clusters)]
    else:
        colors = colors[:len(clusters)]
    col_df = pd.DataFrame(dict(clusters=clusters,
                               r=[i[0] for i in colors],
                               g=[i[1] for i in colors],
                               b=[i[2] for i in colors]))
    data = pd.merge(data, col_df, left_on=feature, right_on='clusters', how='left')

    pngs = {}
    mw, mh = 0, 0

    for batch in slices:
        sub_data = data[data['Batch'] == batch].copy()
        sub_data['cell'] = obj.obs[obj.obs['Batch'] == batch]['cell_id'].astype(int).values
        seg_df = pd.DataFrame(obj.uns[batch]['seg_cell'].flatten(), columns=['cell'])
        seg_df = pd.merge(seg_df, sub_data, how='left', on='cell')
        seg_df.loc[seg_df['r'].isnull(), ['r', 'g', 'b']] = bg_color
        d1, d2 = obj.uns[batch]['seg_cell'].shape
        seg_display = seg_df[['r', 'g', 'b']].to_numpy().reshape(d1, d2, 3).astype(np.uint8)
        seg_display = black_line(obj.uns[batch]['seg_cell'], seg_display, sep=sep, col=line_color)
        w, h = seg_display.shape[:2]
        if compress_factor:
            seg_display = cv2.resize(seg_display, dsize=(int(h * compress_factor), int(w * compress_factor)))
        if batch in angle_dict:
            seg_display = rotate_bound_white_bg(seg_display, angle_dict[batch], bg=bg_color)
        w, h = seg_display.shape[:2]
        if 'FP' in batch:
            seg_display = cv2.resize(seg_display, dsize=(int(h * 500 / 715), int(w * 500 / 715)))
        seg_display = cut_mtx(seg_display, fig_sep)
        w, h = seg_display.shape[:2]
        pngs[batch] = seg_display
        if w > mw:
            mw = w
        if h > mh:
            mh = h
    arr = np.zeros([mw * nrow, mh * ncol, 3], dtype=np.uint8)
    arr[:, :] = bg_color
    n = 0
    aw, ah = arr.shape[:2]
    for batch in pngs:
        seg_display = pngs[batch]
        w, h = seg_display.shape[:2]
        w_ = n // ncol * mw + int((mw - w) / nrow)
        h_ = n % ncol * mh + int((mh - h) / nrow)
        arr[w_:w_ + w, h_:h_ + h, :] = seg_display
        n += 1

    if raw:
        if fname:
            if fname[-4:] == '.svg':
                save_svg(fname, arr, ah, aw)
            else:
                plt.imsave(fname, arr)
        if show:
            plt.figure(figsize=(16, 16))
            plt.imshow(arr)
    else:
        fig = plt.figure(figsize=(16, 17), dpi=dpi)
        gs = gridspec.GridSpec(16, 17)
        if len(colors) < 16:
            ax1 = fig.add_subplot(gs[0, :len(colors)])
        else:
            ax1 = fig.add_subplot(gs[0, :])
        ax1.axis('off')
        ax1.imshow(np.array(colors).astype(np.uint8).reshape(1, len(colors), 3), aspect='auto')
        for ind, cluster in enumerate(clusters):
            ax1.text(ind - 0.5, -0.6, cluster.replace(' ', '\n'), fontdict={'size': str(legend_size)}, ha='left')
        ax2 = fig.add_subplot(gs[1:, :])
        ax2.imshow(arr)
        if not scale:
            ax2.axis('off')
        if title is None:
            title = feature
        ax2.set_title(title)
        if fname is not None:
            plt.savefig(fname)
        if not show:
            plt.close()
    return None


def featureplot_slices_continuous(obj,
                                  feature,
                                  fname=None,
                                  cmap='viridis',
                                  bg_color=(0, 0, 0, 1),
                                  line_color=(0, 0, 0, 1),
                                  title=None,
                                  scale=True,
                                  show=True,
                                  max_=None,
                                  min_=None,
                                  dpi=600,
                                  sep=2,
                                  fig_sep=0,
                                  nrow=1,
                                  ncol=1,
                                  slices=None,
                                  angle_dict=[],
                                  compress_factor=None,
                                  raw=True
                                  ):
    if not slices is None:
        obj = obj[obj.obs['Batch'].isin(slices), :].copy()
    else:
        obj = obj.copy()
        slices = obj.obs['Batch'].unique().tolist()
    cmap = plt.cm.get_cmap(cmap)
    if feature in obj.obs_keys():
        data = obj.obs.loc[:, [feature]]
    else:
        arr = obj[:, feature].X.toarray()
        arr = arr.reshape(arr.shape[0])
        data = pd.DataFrame(arr, index=obj.obs_names, columns=[feature])
    data['Batch'] = obj.obs['Batch']
    data['cell'] = obj.obs['cell_id'].values
    if max_ is not None:
        data.loc[data[feature] > max_, feature] = max_
    else:
        max_ = data[feature].max()
    if min_ is not None:
        data.loc[data[feature] < min_, feature] = min_
    else:
        min_ = data[feature].min()
    data[feature] = (data[feature] - min_) / max_
    data[feature] = data[feature].fillna(0)
    data[['r', 'g', 'b', 'a']] = cmap(data[feature])

    pngs = {}
    mw, mh = 0, 0
    for batch in slices:
        sub_data = data[data['Batch'] == batch].copy()
        sub_data['cell'] = obj.obs[obj.obs['Batch'] == batch]['cell_id'].astype(int).values
        seg_df = pd.DataFrame(obj.uns[batch]['seg_cell'].flatten(), columns=['cell'])

        seg_df = pd.merge(seg_df, sub_data, how='left', on='cell')
        seg_df.loc[seg_df['r'].isnull(), ['r', 'g', 'b', 'a']] = bg_color
        d1, d2 = obj.uns[batch]['seg_cell'].shape
        seg_display = seg_df[['r', 'g', 'b', 'a']].to_numpy().reshape(d1, d2, 4)
        seg_display = black_line(obj.uns[batch]['seg_cell'], seg_display, sep=sep, col=line_color)
        w, h = seg_display.shape[:2]
        if compress_factor:
            seg_display = cv2.resize(seg_display, dsize=(int(h * compress_factor), int(w * compress_factor)))
        if batch in angle_dict:
            seg_display = rotate_bound_white_bg(seg_display, angle_dict[batch], bg=bg_color)
        w, h = seg_display.shape[:2]
        if 'FP' in batch:
            seg_display = cv2.resize(seg_display, dsize=(int(h * 500 / 715), int(w * 500 / 715)))
        seg_display = cut_mtx(seg_display, fig_sep, threshold=1)
        w, h = seg_display.shape[:2]
        pngs[batch] = seg_display
        if w > mw:
            mw = w
        if h > mh:
            mh = h

    arr = np.zeros([mw * nrow, mh * ncol, 4])
    arr[:, :] = bg_color
    n = 0
    aw, ah = arr.shape[:2]
    for batch in pngs:
        seg_display = pngs[batch]
        w, h = seg_display.shape[:2]
        w_ = n // ncol * mw + int((mw - w) / nrow)
        h_ = n % ncol * mh + int((mh - h) / nrow)
        arr[w_:w_ + w, h_:h_ + h, :] = seg_display
        n += 1

    if raw:
        if fname:
            if fname[-4:] == '.svg':
                save_svg(fname, arr, ah, aw)
            else:
                plt.imsave(fname, arr)
        if show:
            plt.figure(figsize=(16, 16))
            plt.imshow(arr)
    else:
        fig = plt.figure(figsize=(16, 17), dpi=dpi)
        gs = gridspec.GridSpec(16, 17)
        ax1 = fig.add_subplot(gs[0, :5])
        ax1.axis('off')
        ax1.imshow(cmap(np.linspace(0, 1, 100)).reshape(1, 100, 4), aspect='auto')
        ax1.text(-0.5, -0.6, f'{min_:.3f}', fontdict={'size': '12'}, ha='left')
        ax1.text(100, -0.6, f'{max_:.3f}', fontdict={'size': '12'}, ha='right')
        ax2 = fig.add_subplot(gs[1:, :])
        ax2.imshow(arr)
        if not scale:
            ax2.axis('off')
        if title is None:
            title = feature
        ax2.set_title(title)
        if fname is not None:
            plt.savefig(fname)
        if not show:
            plt.close()

    return max_, min_


def featureplot_discrete(obj,
                         feature,
                         order=None,
                         colors=None,
                         bg_color=(0, 0, 0),
                         line_color=(0, 0, 0),
                         blank_color=(192, 192, 192),
                         sep=2,
                         slice=None
                         ):
    if colors is None:
        colors = colorlist
    elif isinstance(colors[0], str):
        colors = [Hex_to_RGB(i) for i in colors]
    data = obj.obs[feature].astype(str)
    clusters = data.unique().tolist()
    if order:
        if len(order) < len(clusters):
            clusters = order
            colors = colors[:len(clusters)]
            data[~data.isin(order)] = ''
            clusters.append('')
            colors.append(blank_color)
        else:
            clusters = order
            colors = colors[:len(clusters)]
            colors = colors[:len(clusters)]
    else:
        colors = colors[:len(clusters)]
    col_df = pd.DataFrame(dict(clusters=clusters,
                               r=[i[0] for i in colors],
                               g=[i[1] for i in colors],
                               b=[i[2] for i in colors]))
    data = pd.merge(data, col_df, left_on=feature, right_on='clusters', how='left')
    data['cell'] = obj.obs['cell_id'].astype(int).values
    seg_df = pd.DataFrame(obj.uns[slice]['seg_cell'].flatten(), columns=['cell'])
    seg_df = pd.merge(seg_df, data, how='left', on='cell')
    seg_df.loc[seg_df['r'].isnull(), ['r', 'g', 'b']] = bg_color
    d1, d2 = obj.uns[slice]['seg_cell'].shape
    seg_display = seg_df[['r', 'g', 'b']].to_numpy().reshape(d1, d2, 3).astype(np.uint8)
    seg_display = black_line(obj.uns[slice]['seg_cell'], seg_display, sep=sep, col=line_color)
    return seg_display


def featureplot_single_discrete(obj,
                                feature,
                                fname=None,
                                order=None,
                                colors=None,
                                bg_color=(0,0,0),
                                line_color=(0, 0, 0),
                                blank_color=(192, 192, 192),
                                sep=2,
                                show=True,
                                scale=False,
                                title=None,
                                dpi=600,
                                fig_sep=0,
                                nrow=1,
                                ncol=1,
                                slice=None,
                                angle_dict=None,
                                compress=True,
                                raw=True,
                                legend_size=12,
                                ):
    if angle_dict is None:
        angle_dict = {}
    obj = obj[obj.obs['Batch'] == slice, :].copy()
    if colors is None:
        colors = colorlist
    elif isinstance(colors[0], str):
        colors = [Hex_to_RGB(i) for i in colors]
    data = obj.obs[[feature, 'Batch']]
    data[feature] = data[feature].astype(str)
    clusters = data[feature].unique().tolist()
    if order:
        new_order = []
        new_colors = []
        for i, j in zip(order, colors):
            if i in clusters:
                new_order.append(i)
                new_colors.append(j)
        order = new_order
        colors = new_colors
        if len(order) < len(clusters):
            clusters = order
            colors = colors[:len(clusters)]
            data.loc[~data[feature].isin(order), feature] = ''
            clusters.append('')
            colors.append(blank_color)
        else:
            clusters = order
            colors = colors[:len(clusters)]
            colors = colors[:len(clusters)]
    else:
        colors = colors[:len(clusters)]

    pngs = {}
    mw, mh = 0, 0
    for cluster, color in zip(clusters, colors):
        seg_display = featureplot_discrete(obj,
                                           feature,
                                           order=[cluster],
                                           colors=[color],
                                           slice=slice,
                                           bg_color=bg_color,
                                           line_color=line_color,
                                           sep=sep,
                                           blank_color=blank_color
                                           )
        if slice in angle_dict:
            seg_display = rotate_bound_white_bg(seg_display, angle_dict[slice], bg=bg_color)
        w, h = seg_display.shape[:2]
        if 'FP' in slice:
            seg_display = cv2.resize(seg_display, dsize=(int(h * 500 / 715), int(w * 500 / 715)))
        seg_display = cut_mtx(seg_display, fig_sep)
        w, h = seg_display.shape[:2]
        pngs[cluster] = seg_display
        if w > mw:
            mw = w
        if h > mh:
            mh = h
    arr = np.zeros([mw * nrow, mh * ncol, 3], dtype=np.uint8)
    arr[:, :] = bg_color
    n = 0
    aw, ah = arr.shape[:2]
    for batch in pngs:
        seg_display = pngs[batch]
        w, h = seg_display.shape[:2]
        w_ = n // ncol * mw + int((mw - w) / nrow)
        h_ = n % ncol * mh + int((mh - h) / nrow)
        arr[w_:w_ + w, h_:h_ + h, :] = seg_display
        n += 1

    if raw:
        if fname:
            if fname[-4:] == '.svg':
                save_svg(fname, arr, ah, aw)
            else:
                plt.imsave(fname, arr)
        if show:
            plt.figure(figsize=(16, 16))
            plt.imshow(arr)
    else:
        fig = plt.figure(figsize=(16, 17), dpi=dpi)
        gs = gridspec.GridSpec(16, 17)
        if len(colors) < 16:
            ax1 = fig.add_subplot(gs[0, :len(colors)])
        else:
            ax1 = fig.add_subplot(gs[0, :])
        ax1.axis('off')
        ax1.imshow(np.array(colors).astype(np.uint8).reshape(1, len(colors), 3), aspect='auto')
        for ind, cluster in enumerate(clusters):
            ax1.text(ind - 0.5, -0.6, cluster.replace(' ', '\n'), fontdict={'size': str(legend_size)}, ha='left')
        ax2 = fig.add_subplot(gs[1:, :])
        ax2.imshow(arr)
        if not scale:
            ax2.axis('off')
        if title is None:
            title = feature
        ax2.set_title(title)
        if fname is not None:
            plt.savefig(fname)
        if not show:
            plt.close()
    return None


def get_cell_center(obj):
    obj = obj.copy()
    coor = pd.DataFrame(obj.uns['seg_cell'].flatten(), columns=['cell'])
    coor = coor[coor['cell'] > 0]
    coor[['x', 'y']] = np.argwhere(obj.uns['seg_cell'] > 0)
    obj.obsm['cell_center'] = coor.groupby('cell').mean().values
    return obj


def lr(
        adata,
        lr_pairs: str,
        distance: float,
        name: str = 'cci',
        verbose: bool = True,
):
    """Calculate the proportion of known ligand-receptor co-expression among the neighbouring spots or within spots
    Parameters
    ----------
    adata: AnnData          The data object to scan
    use_lr: str             object to keep the result (default: adata.uns['cci_lr'])
    distance: float         Distance to determine the neighbours (default: closest), distance=0 means within spot
    Returns
    -------
    adata: AnnData          The data object including the results
    """

    df = adata.to_df()

    # expand the LR pairs list by swapping ligand-receptor positions
    # lr_pairs = adata.uns["lr"].copy()
    # lr_pairs += [item.split("_")[1] + "_" + item.split("_")[0] for item in lr_pairs]
    lr_pairs = [lr_pairs]

    # get neighbour spots for each spot according to the specified distance
    coor = adata.obsm['cell_center']
    point_tree = cKDTree(coor)
    neighbours = []

    if distance == 0:
        neighbours = adata.obs_names.to_list()
    else:
        for spot, spot_coor in zip(adata.obs_names, coor):
            n_index = point_tree.query_ball_point(
                spot_coor,
                distance,
            )
            neighbours.append(
                [item for item in df.index[n_index] if not (item == spot)]
            )
    # filter out those LR pairs that do not exist in the dataset
    lr1 = [item.split("_")[0] for item in lr_pairs]
    lr2 = [item.split("_")[1] for item in lr_pairs]
    avail = [
        i for i, x in enumerate(lr1) if lr1[i] in df.columns and lr2[i] in df.columns
    ]
    spot_lr1 = df[[lr1[i] for i in avail]]
    spot_lr2 = df[[lr2[i] for i in avail]]
    if verbose:
        print("Altogether " + str(len(avail)) + " valid L-R pairs")

    # function to calculate mean of lr2 expression between neighbours or within spot (distance==0) for each spot
    def mean_lr2(x):
        # get lr2 expressions from the neighbour(s)
        nbs = spot_lr2.loc[neighbours[df.index.tolist().index(x.name)], :]
        if nbs.shape[0] > 0:  # if neighbour exists
            return (nbs > 0).sum() / nbs.shape[0]
        else:
            return 0

    # mean of lr2 expressions from neighbours of each spot
    nb_lr2 = spot_lr2.apply(mean_lr2, axis=1)

    # check whether neighbours exist
    try:
        nb_lr2.shape[1]
    except:
        raise ValueError("No neighbours found within given distance.")

    # keep value of nb_lr2 only when lr1 is also expressed on the spots
    spot_lr = pd.DataFrame(
        spot_lr1.values * (nb_lr2.values > 0) + (spot_lr1.values > 0) * nb_lr2.values,
        index=df.index,
        columns=[lr_pairs[i] for i in avail],
    ).sum(axis=1)
    adata.obs[name] = spot_lr.values / 2
