import numpy as np
import pandas as pd
import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = pow(2,40).__str__()
import cv2
import matplotlib.pyplot as plt
from tqdm import tqdm
from PIL import Image, ImageDraw
import pywt
from statsmodels.stats.outliers_influence import OLSInfluence
from skimage import filters
import statsmodels.api as sm


class Registration:

    def __init__(self):
        self.ssDNA_path = None
        self.mRNA_path = None
        self.raw_ssDNA = None
        self.raw_mRNA = None

        self.cr_ssDNA_scale = 50
        self.cr_mRNA_scale = 20

        self.ssDNA_e = None
        self.mRNA_e = None
        self.cr_angle = None
        self.cr_factor = None

        self.pr_angle = None
        self.pr_ssDNA = None
        self.ssDNA_track_ys = self.ssDNA_track_xs = None
        self.mRNA_track_ys = self.mRNA_track_xs = None
        self.cr_center = self.adj_center = None
        self.model_xs = self.model_ys = None

        self.matched_ssDNA = self.matched_ssDNA_adj = None

    def __repr__(self):
        t = f"ssDNA Image Registration Object\n" \
            f"ssDNA Image Path: {self.ssDNA_path}\n" \
            f"GEM   Data  Path: {self.mRNA_path}"
        return t

    def load(self,
             ssDNA_path,
             mRNA_path,
             ssDNA_flip=True,
             mRNA_gray_factor=1,
             signal_pbar=None
             ):
        self.ssDNA_path = ssDNA_path
        self.mRNA_path = mRNA_path

        self.raw_ssDNA = cv2.imread(self.ssDNA_path, cv2.IMREAD_GRAYSCALE)
        if ssDNA_flip:
            self.raw_ssDNA = np.flipud(self.raw_ssDNA)
        if signal_pbar:
            signal_pbar.emit(50)
        if mRNA_path[-3:] not in ['png', 'jpg', 'tif']:
            self.raw_mRNA = mRNA_image_load(self.mRNA_path, mRNA_gray_factor)
        else:
            self.raw_mRNA = cv2.imread(self.mRNA_path, cv2.IMREAD_GRAYSCALE)*mRNA_gray_factor
            self.raw_mRNA = self.raw_mRNA.astype(np.uint8)

    def coarse_registration(self,
                            cr_ssDNA_scale=50,
                            cr_mRNA_scale=20,
                            cr_ssDNA_Blur=(9, 9),
                            cr_mRNA_Blur=(13, 13),
                            cr_ssDNA_t='auto',
                            cr_mRNA_t='auto',
                            verbose=True,
                            signal_pbar=None
                            ):
        """
        Rough Registration
        """
        # Maximize the margin & Fitting the ellipse
        self.cr_ssDNA_scale = cr_ssDNA_scale
        ssDNA = cv2.resize(self.raw_ssDNA,
                           (self.raw_ssDNA.shape[1] // self.cr_ssDNA_scale,
                            self.raw_ssDNA.shape[0] // self.cr_ssDNA_scale)
                           )
        ssDNAgb = cv2.GaussianBlur(ssDNA, cr_ssDNA_Blur, 0)

        if cr_ssDNA_t == 'auto':
            thresholds = filters.threshold_multiotsu(ssDNAgb, classes=3)
            cr_ssDNA_t = thresholds[1] - 10
            print('auto recommond cr_ssDNA_t:', cr_ssDNA_t)

        _, ssDNA = cv2.threshold(ssDNAgb, cr_ssDNA_t, 255, cv2.THRESH_BINARY)
        ssDNA_contour = get_max_contour(ssDNA)
        self.ssDNA_e = cv2.fitEllipse(ssDNA_contour)

        if signal_pbar:
            signal_pbar.emit(35)

        self.cr_mRNA_scale = cr_mRNA_scale
        mRNA = cv2.resize(self.raw_mRNA,
                          (self.raw_mRNA.shape[1] // self.cr_mRNA_scale,
                           self.raw_mRNA.shape[0] // self.cr_mRNA_scale)
                          )
        mRNAgb = cv2.GaussianBlur(mRNA, cr_mRNA_Blur, 0)

        if cr_mRNA_t == 'auto':
            thresholds = filters.threshold_multiotsu(mRNAgb, classes=3)
            cr_mRNA_t = thresholds[0]
            print('auto recommond cr_mRNA_t:', cr_mRNA_t)

        _, mRNA = cv2.threshold(mRNAgb, cr_mRNA_t, 255, cv2.THRESH_BINARY)
        mRNA_contour = get_max_contour(mRNA)
        self.mRNA_e = cv2.fitEllipse(mRNA_contour)

        if signal_pbar:
            signal_pbar.emit(55)

        # Calculate the angle and factor for rough registration
        factor = ((self.mRNA_e[1][0] * self.mRNA_e[1][1]) / (self.ssDNA_e[1][0] * self.ssDNA_e[1][1])) ** 0.5
        angle_1 = self.ssDNA_e[2] - self.mRNA_e[2]
        angle_2 = self.ssDNA_e[2] - self.mRNA_e[2] - 180
        contour_1 = contour_rotate(ssDNA_contour, self.ssDNA_e[0], angle_1)
        contour_2 = contour_rotate(ssDNA_contour, self.ssDNA_e[0], angle_2)
        scde = cv2.createShapeContextDistanceExtractor()
        dist_1 = scde.computeDistance(mRNA_contour, contour_1 * factor)
        dist_2 = scde.computeDistance(mRNA_contour, contour_2 * factor)
        self.cr_angle, contour = (angle_1, contour_1) if dist_1 < dist_2 else (angle_2, contour_2)
        self.cr_factor = factor / self.cr_ssDNA_scale * self.cr_mRNA_scale

        if signal_pbar:
            signal_pbar.emit(60)

        self.cr_mRNA = cv2.drawContours(cv2.cvtColor(mRNAgb, cv2.COLOR_GRAY2RGB), [mRNA_contour], -1, (255, 0, 0), 3,
                                        lineType=cv2.LINE_AA)
        if signal_pbar:
            signal_pbar.emit(70)

        cr_ssDNA = cv2.drawContours(cv2.cvtColor(ssDNAgb, cv2.COLOR_GRAY2RGB), [ssDNA_contour], -1, (255, 0, 0), 3, lineType=cv2.LINE_AA)
        cr_ssDNA = cv2.resize(cr_ssDNA,
                              (int(cr_ssDNA.shape[1] * factor),
                               int(cr_ssDNA.shape[0] * factor))
                              )
        self.cr_ssDNA = rotate_bound_white_bg(cr_ssDNA, self.cr_angle, bg=(0, 0, 0))

        # Display
        if verbose:
            plt.figure(figsize=(16, 4))
            _ = plt.hist(ssDNAgb.flatten(), 256, [0, 256])
            plt.axvline(cr_ssDNA_t, color='red')
            plt.title("ssDNA Gradation Histogram")
            plt.figure(figsize=(16, 4))
            _ = plt.hist(mRNAgb.flatten(), 256, [0, 256])
            plt.axvline(cr_mRNA_t, color='red')
            plt.title("mRNA Gradation Histogram")
            mRNA_c = mRNA_contour - mRNA_contour.min(0)
            ssDNA_c = (contour * factor).astype(int)
            ssDNA_c = ssDNA_c - ssDNA_c.min(0)[0] + (mRNA_c.max(0)[0][0], 0)
            h = ssDNA_c.max(0)[0][0]
            w = ssDNA_c.max(0)[0][1] if ssDNA_c.max(0)[0][1] > mRNA_c.max(0)[0][1] else mRNA_c.max(0)[0][1]
            img = np.ones((w, h, 3)) * 255
            cv2.drawContours(img, [ssDNA_c], -1, (0, 0, 255), 3, lineType=cv2.LINE_AA)
            cv2.drawContours(img, [mRNA_c], -1, (0, 255, 0), 3, lineType=cv2.LINE_AA)
            plt.figure(figsize=(8, 8))
            plt.imshow(img)
            plt.title("mRNA - ssDNA Coarse Registration")
            print(f"Fitted ssDNA Ellipse: {self.ssDNA_e}")
            print(f"Fitted mRNA  Ellipse: {self.mRNA_e}")
            print(f"Coarse Rotation     Angle: {self.cr_angle}")
            print(f"Coarse Registration Factor: {self.cr_factor}")

    def precise_rotation(self,
                         pr_ssDNA_scale=5,
                         pr_ssDNA_t0='auto',
                         pr_ssDNA_t1='auto',
                         step=0.1,
                         angle_range=5,
                         verbose=True,
                         only_x = False,
                         only_y = False,
                         signal_pbar=None
                         ):
        '''
        Accurate angle
        '''
        if pr_ssDNA_t0 == 'auto':
            pr_ssDNA_t0 = 63
        if pr_ssDNA_t1 == 'auto':
            pr_ssDNA_t1 = 83
        ssDNA = cv2.resize(self.raw_ssDNA,
                           (self.raw_ssDNA.shape[1] // pr_ssDNA_scale,
                            self.raw_ssDNA.shape[0] // pr_ssDNA_scale)
                           )
        ssDNAt = ssDNA.copy()
        ssDNAt[(ssDNA > pr_ssDNA_t1) | (ssDNA < pr_ssDNA_t0)] = 255
        ssDNAt[(ssDNA <= pr_ssDNA_t1) & (ssDNA >= pr_ssDNA_t0)] = 0

        ldwt_x = []
        ldwt_y = []
        lx = []
        ly = []
        l = np.arange(self.cr_angle - angle_range, self.cr_angle + angle_range, step)
        for angle in tqdm(l):
            img = rotate_bound_white_bg(ssDNAt, angle, bg=(0, 0, 0))
            _, x_dwt = pywt.dwt(img.mean(axis=0), 'db1')
            _, y_dwt = pywt.dwt(img.mean(axis=1), 'db1')
            ldwt_x.append(x_dwt)
            ldwt_y.append(y_dwt)
            std_x = np.std(x_dwt[20:-20])
            std_y = np.std(y_dwt[20:-20])
            lx.append(std_x)
            ly.append(std_y)
        mlx = lx.index(max(lx))
        mly = ly.index(max(ly))
        if only_x:
            l_inx = mlx
        elif only_y:
            l_inx = mly
        else:
            l_inx = np.mean([mlx, mly])
        if abs(mlx-mly) > (0.3/step):
            print(f"Waringing! index of max x's std is {mlx}, while y's is {mly}")
        self.pr_angle = self.cr_angle + l_inx * step - angle_range
        self.pr_ssDNA = rotate_bound_white_bg(self.raw_ssDNA, self.pr_angle, bg=(0, 0, 0))

        if signal_pbar:
            signal_pbar.emit(50)

        if verbose:
            plt.figure(figsize=(16, 4))
            _ = plt.hist(ssDNA.flatten(), 256, [0, 256])
            plt.axvline(pr_ssDNA_t0, color='red')
            plt.axvline(pr_ssDNA_t1, color='red')
            plt.title(f"Scaled ssDNA Gradation Histogram [{pr_ssDNA_t0} {pr_ssDNA_t1}]")

            img = rotate_bound_white_bg(ssDNAt, self.pr_angle, bg=(0, 0, 0))
            _, x_dwt = pywt.dwt(img.mean(axis=0), 'db1')
            _, y_dwt = pywt.dwt(img.mean(axis=1), 'db1')
            plt.figure(figsize=(16, 8))
            ax = plt.subplot(221)
            plt.plot(l, lx)
            ax.set_title('x-axis std')
            ax = plt.subplot(222)
            plt.plot(l, ly)
            ax.set_title('y-axis std')
            ax = plt.subplot(223)
            plt.plot(x_dwt)
            ax.set_title('x-axis dwt at precise angle')
            ax = plt.subplot(224)
            plt.plot(y_dwt)
            ax.set_title('y-axis dwt at precise angle')
            plt.figure(figsize=(16, 16))
            plt.imshow(img, "gray")
            plt.title(f"Rotated ssDNA[{pr_ssDNA_t0}: {pr_ssDNA_t1}] Image")
            plt.figure(figsize=(16, 16))
            plt.imshow(self.pr_ssDNA, "gray")
            plt.title("Rotated ssDNA Image")
            print(f"Precise Rotation Angle: {self.pr_angle}")

    def track_detection(self,
                        filter_kernel=None,
                        tf_ssDNA_scale=5,
                        tf_ssDNA_t0='auto',
                        tf_ssDNA_t1='auto',
                        morph_rect_size=350,
                        verbose=True,
                        signal_pbar=None
                        ):
        '''
        Identify track line
        '''
        if tf_ssDNA_t0 == 'auto':
            tf_ssDNA_t0 = 48
        if tf_ssDNA_t1 == 'auto':
            tf_ssDNA_t1 = 62
        self.tf_ssDNA_scale = tf_ssDNA_scale
        if filter_kernel is None:
            filter_kernel = -np.ones((201, 201), dtype=np.int8)
            filter_kernel[100] = 0
            filter_kernel[100, 20:180] = 40
        ssDNA = cv2.resize(self.pr_ssDNA,
                           (self.pr_ssDNA.shape[1] // tf_ssDNA_scale,
                            self.pr_ssDNA.shape[0] // tf_ssDNA_scale)
                           )
        ssDNAt = ssDNA.copy()
        ssDNAt[(ssDNA > tf_ssDNA_t1) | (ssDNA < tf_ssDNA_t0)] = 0
        ssDNAt[(ssDNA <= tf_ssDNA_t1) & (ssDNA >= tf_ssDNA_t0)] = 255

        ssDNAout_y = cv2.filter2D(ssDNAt, -1, filter_kernel)
        kernel_y = cv2.getStructuringElement(cv2.MORPH_RECT, (morph_rect_size, 1))
        ssDNAex_y = cv2.morphologyEx(ssDNAout_y, cv2.MORPH_OPEN, kernel_y)
        self.ssDNA_track_ys = filter_n(np.where(ssDNAex_y.mean(axis=1))[0]) * tf_ssDNA_scale

        ssDNAout_x = cv2.filter2D(ssDNAt, -1, filter_kernel.T)
        kernel_x = cv2.getStructuringElement(cv2.MORPH_RECT, (1, morph_rect_size))
        ssDNAex_x = cv2.morphologyEx(ssDNAout_x, cv2.MORPH_OPEN, kernel_x)
        self.ssDNA_track_xs = filter_n(np.where(ssDNAex_x.mean(axis=0))[0]) * tf_ssDNA_scale

        self.mRNA_track_ys = get_mRNA_track_lines(self.raw_mRNA.mean(axis=1))
        self.mRNA_track_xs = get_mRNA_track_lines(self.raw_mRNA.mean(axis=0))

        # Transform the center
        ssDNA_old_center = (self.ssDNA_e[0][0] * self.cr_ssDNA_scale, self.ssDNA_e[0][1] * self.cr_ssDNA_scale)
        self.cr_center = rotate_point(self.raw_ssDNA, self.pr_angle, ssDNA_old_center)

        if signal_pbar:
            signal_pbar.emit(70)

        ssDNA = cv2.cvtColor(cv2.resize(self.pr_ssDNA,
                                        (self.pr_ssDNA.shape[1] // tf_ssDNA_scale,
                                         self.pr_ssDNA.shape[0] // tf_ssDNA_scale)
                                        ), cv2.COLOR_GRAY2RGB)
        self.pr_ssDNA_track = Image.fromarray(ssDNA)
        img_draw = ImageDraw.Draw(self.pr_ssDNA_track)
        fill = (0, 255, 0)
        for x in self.ssDNA_track_xs:
            x = x // tf_ssDNA_scale
            img_draw.line((x, 0) + (x, ssDNA.shape[0]), fill=fill, width=3)
        for y in self.ssDNA_track_ys:
            y = y // tf_ssDNA_scale
            img_draw.line((0, y) + (ssDNA.shape[1], y), fill=fill, width=3)

        rgb = cv2.cvtColor(self.raw_mRNA, cv2.COLOR_GRAY2RGB)
        self.mRNA_track = Image.fromarray(rgb)
        img_draw = ImageDraw.Draw(self.mRNA_track)
        fill2 = (255, 0, 0)
        r = 8
        for x in self.mRNA_track_xs:
            for y in self.mRNA_track_ys:
                img_draw.ellipse((x - r, y - r, x + r, y + r), fill=fill2)
        xs = (self.ssDNA_track_xs - self.cr_center[0]) * self.cr_factor + \
             (self.mRNA_e[0][0] * self.cr_mRNA_scale)
        ys = (self.ssDNA_track_ys - self.cr_center[1]) * self.cr_factor + \
             (self.mRNA_e[0][1] * self.cr_mRNA_scale)
        for x in xs:
            img_draw.line((x, 0) + (x, self.raw_mRNA.shape[0]), fill=fill, width=5)
        for y in ys:
            img_draw.line((0, y) + (self.raw_mRNA.shape[1], y), fill=fill, width=5)

        if verbose:
            plt.figure(figsize=(16, 4))
            _ = plt.hist(ssDNA.flatten(), 256, [0, 256])
            plt.axvline(tf_ssDNA_t0, color='red')
            plt.axvline(tf_ssDNA_t1, color='red')
            plt.title(f"Scaled ssDNA Gradation Histogram [{tf_ssDNA_t0} {tf_ssDNA_t1}]")

            plt.figure(figsize=(16, 16))
            plt.imshow(ssDNAt, "gray")
            plt.title(f"Rotated ssDNA[{tf_ssDNA_t0}: {tf_ssDNA_t1}] Image")
            plt.figure(figsize=(16, 8))
            ax = plt.subplot(221)
            plt.imshow(ssDNAout_y, "gray")
            ax.set_title('Horizontal-Convolution')
            ax = plt.subplot(222)
            plt.imshow(ssDNAout_x, "gray")
            ax.set_title('Vertical-Convolution')
            ax = plt.subplot(223)
            plt.imshow(ssDNAex_y, "gray")
            ax.set_title('Horizontal-morphologyEx')
            ax = plt.subplot(224)
            plt.imshow(ssDNAex_x, "gray")
            ax.set_title('Vertical-morphologyEx')

            plt.figure(figsize=(16, 16))
            plt.imshow(self.pr_ssDNA_track)
            plt.title("Rotated ssDNA with Track Lines")

            plt.figure(figsize=(16, 16))
            plt.imshow(self.mRNA_track)
            plt.title("mRNA with Track Points")

    def precise_registration(self,
                             adj_x=0,
                             adj_y=0,
                             rsq=0.99999,
                             verbose=True,
                             signal_pbar=None
                             ):
        '''
        Accurate Registration
        :param adj_x:
        :param adj_y:
        :param rsq:
        :param verbose:
        :return:
        '''
        self.adj_center = self.cr_center[0] + adj_x, self.cr_center[1] + adj_y
        xs = (self.ssDNA_track_xs - self.adj_center[0]) * self.cr_factor + \
             (self.mRNA_e[0][0] * self.cr_mRNA_scale)
        ys = (self.ssDNA_track_ys - self.adj_center[1]) * self.cr_factor + \
             (self.mRNA_e[0][1] * self.cr_mRNA_scale)
        xs = xs[xs >= -100]
        xs = xs[xs < self.raw_mRNA.shape[0]]
        ys = ys[ys >= -100]
        ys = ys[ys < self.raw_mRNA.shape[1]]

        if verbose:
            rgb = cv2.cvtColor(self.raw_mRNA, cv2.COLOR_GRAY2RGB)
            img_pil = Image.fromarray(rgb)
            img_draw = ImageDraw.Draw(img_pil)
            fill = (0, 255, 0)
            fill2 = (255, 0, 0)
            r = 8
            for x in self.mRNA_track_xs:
                for y in self.mRNA_track_ys:
                    img_draw.ellipse((x - r, y - r, x + r, y + r), fill=fill2)
            for x in xs:
                img_draw.line((x, 0) + (x, self.raw_mRNA.shape[0]), fill=fill, width=5)
            for y in ys:
                img_draw.line((0, y) + (self.raw_mRNA.shape[1], y), fill=fill, width=5)
            plt.figure(figsize=(16, 16))
            plt.imshow(img_pil)
            plt.title("mRNA Before Matching")

        if signal_pbar:
            signal_pbar.emit(40)

        neighbor_xs = get_neighbor(xs, self.mRNA_track_xs)
        neighbor_ys = get_neighbor(ys, self.mRNA_track_ys)
        matched_xs, matched_neighbor_xs, self.model_xs = match_track_lines(xs, neighbor_xs, rsq)
        matched_ys, matched_neighbor_ys, self.model_ys = match_track_lines(ys, neighbor_ys, rsq)

        adj_xs = matched_xs * self.model_xs.params[1] + self.model_xs.params[0]
        adj_ys = matched_ys * self.model_ys.params[1] + self.model_ys.params[0]

        rgb = cv2.cvtColor(self.raw_mRNA, cv2.COLOR_GRAY2RGB)
        self.matched_mRNA = Image.fromarray(rgb)
        img_draw = ImageDraw.Draw(self.matched_mRNA)
        fill = (0, 255, 0)
        fill2 = (255, 0, 0)
        r = 8
        for x in self.mRNA_track_xs:
            for y in self.mRNA_track_ys:
                img_draw.ellipse((x - r, y - r, x + r, y + r), fill=fill2)
        for x in adj_xs:
            img_draw.line((x, 0) + (x, self.raw_mRNA.shape[0]), fill=fill, width=5)
        for y in adj_ys:
            img_draw.line((0, y) + (self.raw_mRNA.shape[1], y), fill=fill, width=5)

        if signal_pbar:
            signal_pbar.emit(70)

        if verbose:
            plt.figure(figsize=(16, 16))
            plt.imshow(self.matched_mRNA)
            plt.title("mRNA After Matching")

        x0 = (0 - self.model_xs.params[0]) / self.model_xs.params[1]
        x0 = (x0 - (self.mRNA_e[0][0] * self.cr_mRNA_scale)) / self.cr_factor + self.adj_center[0]
        x1 = (self.raw_mRNA.shape[1] - self.model_xs.params[0]) / self.model_xs.params[1]
        x1 = (x1 - (self.mRNA_e[0][0] * self.cr_mRNA_scale)) / self.cr_factor + self.adj_center[0]

        y0 = (0 - self.model_ys.params[0]) / self.model_ys.params[1]
        y0 = (y0 - (self.mRNA_e[0][1] * self.cr_mRNA_scale)) / self.cr_factor + self.adj_center[1]
        y1 = (self.raw_mRNA.shape[0] - self.model_ys.params[0]) / self.model_ys.params[1]
        y1 = (y1 - (self.mRNA_e[0][1] * self.cr_mRNA_scale)) / self.cr_factor + self.adj_center[1]

        matched_ssDNA = self.pr_ssDNA[int(y0):int(y1) + 1, int(x0):int(x1) + 1]
        self.matched_ssDNA = cv2.resize(matched_ssDNA, (self.raw_mRNA.shape[1], self.raw_mRNA.shape[0]))

        if verbose:
            plt.figure(figsize=(16, 16))
            plt.imshow(matched_ssDNA)
            plt.title("ssDNA After Matching")
            plt.figure(figsize=(16, 16))
            plt.imshow(self.matched_ssDNA)
            plt.title("ssDNA After Matching&Scaling")


def get_neighbor(xs, ys):
    new_ys = []
    for i in xs:
        ind = np.argmax(-np.abs((ys - i)))
        new_ys.append(ys[ind])
    return np.array(new_ys)


def match_track_lines(xs, ys, rsq=0.999):
    xs = xs.copy()
    ys = ys.copy()
    Xs = sm.add_constant(xs)
    model = sm.OLS(ys, Xs).fit()
    print('init      r2:', model.rsquared)
    while model.rsquared < rsq:
        cooks_d = OLSInfluence(model).cooks_distance[0]
        ind = np.argmax(cooks_d)
        xs = np.delete(xs, ind)
        ys = np.delete(ys, ind)
        Xs = sm.add_constant(xs)
        model = sm.OLS(ys, Xs).fit()
        print('iteration r2:', model.rsquared)
    print('final     r2:', model.rsquared, '\n')
    return xs, ys, model


def get_max_contour(img):
    contours, hierarchy = cv2.findContours(img, cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)
    max_contour = contours[0]
    for i in contours[1:]:
        if i.shape[0] > max_contour.shape[0]:
            max_contour = i
    return max_contour


def contour_rotate(contour, center, angle):
    angle = angle / 180 * np.pi
    rot_mat = np.array([[np.cos(angle), -np.sin(angle)],
                        [np.sin(angle), +np.cos(angle)]])
    contour = np.dot((contour - center), rot_mat) + center
    contour = contour.astype(np.int32)
    return contour


def rotate_point(image, angle, point):
    (h, w) = image.shape[:2]
    M, nW, nH = get_rotation_matrix(h, w, angle)
    point = np.array([point[0], point[1], 1.0])
    new_point = np.dot(M, point)
    return new_point


def get_rotation_matrix(h, w, angle):
    (cX, cY) = (w // 2, h // 2)
    M = cv2.getRotationMatrix2D((cX, cY), angle, 1.0)
    cos = np.abs(M[0, 0])
    sin = np.abs(M[0, 1])
    nW = int((h * sin) + (w * cos))
    nH = int((h * cos) + (w * sin))
    M[0, 2] += (nW / 2) - cX
    M[1, 2] += (nH / 2) - cY
    return M, nW, nH


def rotate_bound_white_bg(image, angle, bg=(0, 0, 0)):
    (h, w) = image.shape[:2]
    M, nW, nH = get_rotation_matrix(h, w, angle)
    img = cv2.warpAffine(image, M, (nW, nH), borderValue=bg)
    return img


def filter_n(l, n=1):
    out = []
    last = -999
    flag = False
    for i in l:
        if i - last > n:
            out.append(i)
            flag = False
        elif not flag:
            out = out[:-1]
            flag = True
        last = i
    return np.array(out)


def get_mRNA_track_lines(l, width=3):
    lines = []
    cont = 0
    for ind, i in enumerate(l):
        if i == 0:
            cont += 1
        elif i > 0:
            if cont == width:
                lines.append(ind - (width / 2 + 0.5))
            cont = 0
    return np.array(lines)


def mRNA_image_load(path, gray_factor):
    cbs = pd.read_csv(path, sep='\t',comment = "#") # skiprows = 6
    cbs['coor'] = cbs['x'].astype(str) + '_' + cbs['y'].astype(str)
    #cbs2 = cbs[['coor', 'UMICount']] if 'UMICount' in cbs.columns elif cbs[['coor', 'MIDCount']] else cbs[['coor', 'MIDCounts']]
    if 'UMICount' in cbs.columns:
        cbs2 = cbs[['coor', 'UMICount']]
    elif 'MIDCounts' in cbs.columns:
        cbs2 = cbs[['coor', 'MIDCounts']]
    elif 'MIDCount' in cbs.columns:
        cbs2 = cbs[['coor', 'MIDCount']]

    cbs2 = cbs2['coor'].value_counts().reset_index()
    cbs2['x'] = cbs2['index'].str.split('_', expand=True)[0].astype(int)
    cbs2['y'] = cbs2['index'].str.split('_', expand=True)[1].astype(int)
    cbs2['x'] -= cbs2['x'].min()
    cbs2['y'] -= cbs2['y'].min()
    h = cbs2['x'].max()
    w = cbs2['y'].max()
    mtx = cbs2.pivot_table(values='coor', index='x', columns='y', aggfunc='sum')
    for i in range(h):
        if i not in mtx.index:
            mtx.loc[i] = 0
    for j in range(w):
        if j not in mtx.columns:
            mtx.loc[:, j] = 0
    mtx = mtx.fillna(0).sort_index().sort_index(axis=1)
    mtx = (mtx * gray_factor).astype(np.uint8).values
    return mtx
