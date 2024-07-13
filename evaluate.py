import wfdb
import numpy as np
import os

def read_annot(idx):
    mapping = {
        1: 101,
        2: 106,
        3: 115,
        4: 207,
        5: 233,
        6: 102,
        7: 108,
        8: 119,
        9: 209,
        10: 234
    }
    annot_name = mapping[idx]
    filepath = os.path.join("reference", f"{annot_name}")
    start_pos = 0
    fs = 360
    # 读取注释文件
    annotation = wfdb.rdann(filepath, 'atr', sampfrom=start_pos*fs, sampto=fs*(60+start_pos), shift_samps=True)

    V_indices = [i for i, symbol in enumerate(annotation.symbol) if symbol == 'V']
    V_samples = annotation.sample[V_indices]

    A_indices = [i for i, symbol in enumerate(annotation.symbol) if symbol == 'A']
    A_samples = annotation.sample[A_indices]
    
    return annotation.sample/fs, V_samples/fs, A_samples/fs

def compare_annot(truth, predict, tolerance=0.1):
    true_positives = 0
    false_positives = 0
    false_negatives = 0

    for t in truth:
        # 查找在tolerance范围内的预测心跳
        matching_predicts = np.abs(predict - t) < tolerance
        if np.any(matching_predicts):
            true_positives += 1
            predict = predict[~matching_predicts]
        else:
            false_negatives += 1

    false_positives = len(predict)

    return true_positives, false_positives, false_negatives
