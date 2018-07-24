import pickle
from sklearn.externals import joblib
from sklearn import datasets
from skimage.feature import hog
from sklearn.svm import LinearSVC
import numpy as np
from collections import Counter
import cv2
from skimage.feature import hog
from matplotlib import pyplot as plt
from sklearn.cluster import KMeans
import pandas as pd
from os import listdir
from os.path import isfile, join


def pred_img(image_name):

    def inner_crop(img, pad=10):
        img_h, img_w = src_img.shape
        w = img_w - 2 * pad
        h = img_h - 2 * pad
        return img[pad:pad+h, pad:pad+w]

    # Read the input image
    src_img = cv2.imread('seg/' + image_name + '.jpg', cv2.IMREAD_GRAYSCALE)
    src_img = inner_crop(src_img, 20)

    # Reference para
    img_h, img_w = src_img.shape

    blur_img = cv2.GaussianBlur(src_img, (5, 5), 0)


    # Threshold the image
    _, th_img = cv2.threshold(blur_img, 90, 255, cv2.THRESH_BINARY_INV)

    # Find contours in the image
    _, ctrs, _ = cv2.findContours(th_img.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    # Get rectangles contains each contour
    rects = [cv2.boundingRect(ctr) for ctr in ctrs]
    print(len(rects), ' rectangle found.')

    '''
    Remove rects which have too long width or height.
    '''
    print('#---- Pre-Sift Rects by Size ----#')
    def siftBySize(rect):
        if rect[2] > img_w / 4 or rect[3] > img_h / 4:
            # print('delete rect: ', rect)
            return False
        else:
            return True
        
    rects = list(filter(siftBySize, rects))
    print(len(rects), ' rects last.')
    print('#---- Pre-Sift Rects by Size ----#')

    def save_drawn_rects(img, rects, name='rect'):
        for idx, rect in enumerate(rects):
                cv2.rectangle(img, (rect[0], rect[1]), (rect[0] + rect[2], rect[1] + rect[3]), 130, 3)

        cv2.imwrite('inter/' + name + '.jpg', img)


    print('#---- Line Cluster & Remove Outliers ----#')
    def get_cluster_label(rects):
        estimator = KMeans(n_clusters=9)
        estimator.fit([[rect[1]+rect[3]/2, 0] for rect in rects])
        return estimator.labels_

    labels = get_cluster_label(rects)

    def get_label_of_outlier(labels):
        uni, unic = np.unique(labels, return_counts=True)
        # print("uni: ", uni)
        # print("unic: ", unic)
        # labels_of_outlier = np.append(uni[unic < 8], uni[unic > 100])
        labels_of_outlier = uni[unic < 4]
        # print("labels_of_outlier: ", labels_of_outlier)
        return labels_of_outlier

    labels_of_outlier = get_label_of_outlier(labels)

    def is_outlier(idx):
        label = labels[idx]
        return (labels_of_outlier == label).any()

    loop = 0
    while len(labels_of_outlier) != 0:
        loop += 1
        # print('REMOVE ', loop, ': remove outlier in clusters...')
        rects = [rect for idx, rect in enumerate(rects) if not is_outlier(idx)]
        print('', len(rects), ' rects last.\n--')
        labels = get_cluster_label(rects)
        labels_of_outlier = get_label_of_outlier(labels)

    # print('[Finished.]')
    print(len(rects), ' rects last.')
    print('#---- Line Cluster & Remove Outliers ----#')

    print('#---- Sort Lines ----#')

    uni = np.unique(labels)
    idxes = np.unique(labels, return_index=True)[1]
    lines_in_label = [labels[idx] for idx in sorted(idxes, reverse=True)]
    print('lines_in_label: ', lines_in_label)

    rects_in_line = []
    rects_ = np.asarray(rects)
    for i in range(9):
        print('Sort line', i)
        line = sorted(rects_[labels == lines_in_label[i]], key=lambda item: item[0])
        rects_in_line.append(line)
        
    print('#---- Sort Lines ----#')

    def digit_cluster(rects_in_line, line_num):
        line = rects_in_line[line_num]
        print('--- Clustering line', line_num, '...')
        # print(line)
        n_in_line = [8, 11, 18, 8, 11, 18, 8, 11, 18]

        def get_digit_cluster_label(rects, n):
            estimator = KMeans(n_clusters=n)
            estimator.fit([[rect[0]+rect[2]/2, 0] for rect in rects])
            return estimator.labels_
        
        labels = get_digit_cluster_label(line, n_in_line[line_num])
        # print('labels: ', labels)

        def get_label_of_overlayed(labels):
            uni, unic = np.unique(labels, return_counts=True)
            # print("uni: ", uni)
            # print("unic: ", unic)
            labels_of_overlayed = uni[unic > 1]
            # print("labels_of_overlayed: ", labels_of_overlayed)
            return labels_of_overlayed

        labels_of_overlayed = get_label_of_overlayed(labels)
        # print('labels_of_overlayed: ', labels_of_overlayed)

        def is_overlayed(idx):
            label = labels[idx]
            return (labels_of_overlayed == label).any()

        def merge_rects(rects):
            x = min([rect[0] for rect in rects])
            y = min([rect[1] for rect in rects])
            w = max([rect[0] - x + rect[2] for rect in rects])
            h = max([rect[1] - y + rect[3] for rect in rects])
            return np.asarray([x, y, w, h])
        
        
        print('Merge Rects...')
        '''
        Two Methods of merging
        '''

        # No.1 Naive 1 Round
        print('[Naive 1 Round Mergeing]')
        new_rects = []
        for label in labels_of_overlayed:
            # print('- Merging rect with label', label, '...')
            overlayeds = np.asarray(line)[labels == label]
            # print('Merge:')
            # print(overlayeds)
            # print('into:')
            new_rect = merge_rects(overlayeds)
            # print(new_rect)
            new_rects.append(new_rect)
        line = [rect for idx, rect in enumerate(line) if not is_overlayed(idx)]
        line = line + new_rects
        line.sort(key=lambda item: item[0])


        # No.2 Iteration
        # loop = 0
        # while len(labels_of_overlayed) != 0:
        #     loop += 1
        #     print('### REMOVE ', loop, ': merging digits...')


        #     print('### ', len(line), ' rects last.')
        #     labels = get_cluster_label(rects)
        #     labels_of_outlier = get_label_of_outlier(labels)

        # print('######## Finished.')
        
        # print('Line', line_num, ' clustered.')
        # print(line)
        return line

    print('#---- Digit Cluster ----#')
    lines = []
    n_in_line = [8, 11, 18, 8, 11, 18, 8, 11, 18]
    for line_num in range(9):
        if len(rects_in_line[line_num]) > n_in_line[line_num]:
            new_line = digit_cluster(rects_in_line, line_num)
            save_drawn_rects(src_img.copy(), new_line, 'line' + str(line_num))
            lines.append(new_line)
        else:
            lines.append(rects_in_line[line_num])
    print('#---- Digit Cluster ----#')

    clf = joblib.load("model/digits_cls.pkl")
    res_img = src_img.copy()
    preds = []
    # def pred_line(line, clf):
    for idx, line in enumerate(lines):
        # For each rectangular region, calculate HOG features and predict
        # the digit using Linear SVM.
        line_preds = []
        for rect in line:
            cv2.rectangle(res_img, (rect[0], rect[1]), (rect[0] + rect[2], rect[1] + rect[3]), 130, 3)
            # Make the rec tangular region around the digit
            leng = int(rect[3] * 1.6)
            pt1 = int(rect[1] + rect[3] // 2 - leng // 2)
            pt2 = int(rect[0] + rect[2] // 2 - leng // 2)
            roi = th_img[pt1:pt1+leng, pt2:pt2+leng]
        #     rois.append(roi)
        # for roi in rois:
            if roi.size <= 0: 
                continue
            # Resize the image
            roi = cv2.resize(roi, (28, 28), interpolation=cv2.INTER_AREA)
            roi = cv2.dilate(roi, (3, 3))

            # Calculate the HOG features
            roi_hog_fd = hog(roi, orientations=9, pixels_per_cell=(14, 14), cells_per_block=(1, 1), visualise=False)
            nbr = clf.predict(np.array([roi_hog_fd], 'float64'))
            line_preds.append(nbr[0])
            # print('predict: ', nbr)
            cv2.putText(res_img, str(int(nbr[0])), (rect[0], rect[1]),cv2.FONT_HERSHEY_DUPLEX, 2, 150, 3)

    #     cv2.imwrite('inter/line' + str(idx) + 'pred.jpg', num_img)
        line_preds = line_preds[0:n_in_line[idx]]
        preds.append(''.join(str(e) for e in line_preds))
        
    cv2.imwrite('submit/' + image_name + 'pred.jpg', res_img)
    np.asarray(preds)
    submit = pd.DataFrame(preds)
    submit.to_excel('submit/' + image_name + '.xlsx', header=['number'], index=False)

if __name__ == '__main__':
    mypath = 'seg'
    image_names = [f[:8] for f in listdir(mypath) if isfile(join(mypath, f))]
    # pred_img('15331009')
    for image_name in image_names:
        print('------------------------------------------', image_name)
        pred_img(image_name)