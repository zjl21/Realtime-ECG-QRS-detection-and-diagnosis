import numpy as np
import pandas as pd
from intervaltree import Interval, IntervalTree
from tqdm import tqdm
import matplotlib.pyplot as plt
from detect_R import *
from evaluate import *
import configargparse

parser = configargparse.ArgumentParser(description='Realtime R peak detection and diagnosis based on Pan-Tompkins algorithm')

parser.add_argument('-c', '--config', required=False, is_config_file=True, help='Path to config file.')

# 基本参数
parser.add_argument('--fs', type=int, default=200, help='sampling frequency of the data')
parser.add_argument('--data_idx', type=int, nargs='+', required=True, help='index of the data to be processed')

# 可视化参数
parser.add_argument('--visible', action='store_true', help='whether to visualize the detection process and result')
parser.add_argument('--show_result', action='store_true', help='whether to show the final result')
parser.add_argument('--evaluate', action='store_true', help='whether to evaluate the detection result')
parser.add_argument('--normal_color', type=str, default='green', help='color to display normal heart rate')
parser.add_argument('--tachycardia_color', type=str, default='red', help='color to display tachycardia')
parser.add_argument('--bradycardia_color', type=str, default='blue', help='color to display bradycardia')

# 可变参数
parser.add_argument('--QRS_for_bpm', type=int, default=9, help='number of QRS for bpm calculation')
parser.add_argument('--n_overlap', type=int, default=5, help='overlap data number for each round of real-time peak detection')
parser.add_argument('--chunk_size', type=int, default=100, help='number of data points for each round of real-time peak detection')

# 预处理参数
parser.add_argument('--lc', type=float, default=6.5, help='lower cutoff frequency for high-pass filter for preprocessing')
parser.add_argument('--hc', type=float, default=22, help='higher cutoff frequency for low-pass filter for preprocessing')
parser.add_argument('--order', type=int, default=5, help='order of the filter for proprocessing')
parser.add_argument('--wl', type=int, default=31, help='window length for moving average for preprocessing')

# 阈值迭代参数
parser.add_argument('--param_N', type=float, default=0.75, help='Pan-Tompkins algorithm parameter, see README for details')
parser.add_argument('--param_T', type=float, default=0.5, help='Pan-Tompkins algorithm parameter, see README for details')
parser.add_argument('--thr_T', type=int, default=46, help='Pan-Tompkins algorithm parameter, see README for details')
parser.add_argument('--iter_peak', type=float, default=0.125, help='Pan-Tompkins algorithm parameter, see README for details')
parser.add_argument('--iter_re_peak', type=float, default=0.25, help='Pan-Tompkins algorithm parameter, see README for details')
parser.add_argument('--thr_T_quotient', type=float, default=0.25, help='Pan-Tompkins algorithm parameter, see README for details')

# 早搏诊断参数
parser.add_argument('--thr_width', type=float, default=0.3, help='premature beat diagnosis, see README for details')
parser.add_argument('--thr_quotient', type=float, default=1.5, help='premature beat diagnosis, see README for details')

parser.add_argument('--eval_result_path', type=str, default='result.xlsx', help='the path to save the result for detection evaluation')
parser.add_argument('--result_path', type=str, default='output', help='the folder path to save the result for detection')

args = parser.parse_args()

fs = args.fs
data_idxs = args.data_idx

visible = args.visible
show_result = args.show_result
evaluate = args.evaluate
normal_color = args.normal_color
tachycardia_color = args.tachycardia_color
bradycardia_color = args.bradycardia_color

QRS_for_bpm = args.QRS_for_bpm
n_overlap = args.n_overlap
chunk_size = args.chunk_size
lc = args.lc
hc = args.hc
order = args.order
wl = args.wl
param_N = args.param_N
param_S = 1 - param_N
param_T = args.param_T
iter_peak = args.iter_peak
iter_NS = 1 - iter_peak
iter_re_peak = args.iter_re_peak
iter_re_NS = 1 - iter_re_peak
thr_T = args.thr_T
thr_T_quotient = args.thr_T_quotient
thr_width = args.thr_width
thr_quotient = args.thr_quotient
eval_result_path = args.eval_result_path
result_path = args.result_path



if evaluate:
    tmp_num = len(data_idxs)+2
    accuracy = np.zeros(tmp_num)
    sensitivity = np.zeros(tmp_num)
    premature_accuracy = np.zeros(tmp_num)
    premature_sensitivity = np.zeros(tmp_num)
    f1 = np.zeros(tmp_num)
    f1_premature = np.zeros(tmp_num)
    tp_all = np.zeros(tmp_num)
    fp_all = np.zeros(tmp_num)
    fn_all = np.zeros(tmp_num)
    premature_tp_all = np.zeros(tmp_num)
    premature_fp_all = np.zeros(tmp_num)
    premature_fn_all = np.zeros(tmp_num)
    f1_mean = np.array([])
    f1_outlier_free_mean = np.array([])

if evaluate or show_result:
    os.makedirs(result_path, exist_ok=True)
    # 可视化设置
    colors = {
        "Normal": normal_color,
        "Tachycardia": tachycardia_color,
        "Bradycardia": bradycardia_color,
        "None": "white"
    }


def calculate_heart_rate(QRS,fs=200):
    global Normal, Tachycardia, Bradycardia
    RR_intervals = np.diff(QRS)
    RR_mean = np.mean(RR_intervals)
    bpm = 60 / RR_mean * fs
    if len(QRS) >= 9:
    # 窦性心率诊断
        RR_mean = np.mean(RR_intervals[-8:])
        if  RR_mean < 0.6*fs:
            status = 'Tachycardia'                # 心动过速
            Tachycardia.add(Interval(QRS[0], QRS[-1]))
        elif (RR_mean > 1.2*fs) or RR_intervals[-1] > 1.5*fs:
            status = 'Bradycardia'                # 心动过缓
            Bradycardia.add(Interval(QRS[0], QRS[-1]))
        else:
            status = 'Normal'
            Normal.add(Interval(QRS[0], QRS[-1]))
    else:
        status = 'None'
    return bpm, RR_intervals, status

def diagnose(QRS,width,fs=200,thr_width=thr_width,thr_quotient=thr_quotient):
        if len(QRS) >= 3:
            # 早搏诊断
            RR_intervals = np.diff(QRS)
            quotient = RR_intervals[-1]/RR_intervals[0]
            if (quotient > thr_quotient) or (width[1] > thr_width*fs):
                premature.append(QRS[1])
                premature_new.append(QRS[1])

def calculate_width(QRS):
    width = []
    for i in range(len(QRS)):
        left = QRS[i]
        right = QRS[i]
        if data[QRS[i]] > 0:
            while (left > 2) and ((data[left] > data[left-1]) or (data[left-1] > data[left-2])):
                left -= 1
            while (right < len(data)-2) and ((data[right] > data[right+1]) or (data[right+1] > data[right+2])):
                right += 1
            width = np.append(width, right - left + 1)
        if data[QRS[i]] < 0:
            while (left > 2) and ((data[left] < data[left-1]) or (data[left-1] < data[left-2])):
                left -= 1
            while (right < len(data)-2) and ((data[right] < data[right+1]) or (data[right+1] < data[right+2])):
                right += 1
            width = np.append(width, right - left + 1)
    return width


for data_idx in tqdm(data_idxs):
    # 初始化
    data_size_previous = fs * 8
    data_all = np.array([])
    QRS = np.array([])
    QRS_old = np.array([])
    QRS_new = np.array([])
    QRS_new_all = np.array([])
    Bradycardia = IntervalTree()
    Tachycardia = IntervalTree()
    Normal = IntervalTree()
    Bradycardia_range = np.array([])
    Tachycardia_range = np.array([])
    Normal_range = np.array([])
    width_all = np.array([])
    premature = []
    premature_new = []
    plot_first = True
    peak_first = True

    # 读取数据
    data_all = pd.read_csv(f'data/{data_idx}.txt', header=None)
    data_all = np.squeeze(data_all)


    if visible or show_result:
        # 可视化设定
        fig, ax = plt.subplots(figsize=(40, 6))
        y_min = np.min(data_all) - 0.1
        y_max = np.max(data_all) + 0.2
        ax.set_xlim([0, len(data_all)/fs])
        ax.set_ylim([np.min(data_all), np.max(data_all)])
        ax.set_xlabel('Time(s)')
        ax.set_ylabel('Amplitude')
        ax.set_title(f'R peak detection for data {data_idx}')
        ax.autoscale(enable=False)


    for ii in range(chunk_size, len(data_all), chunk_size):
        data = np.array(data_all[:ii])
        QRS_new_all = []
        if visible:
            ax.set_xlim([0, len(data_all)/fs])
            ax.set_ylim([y_min, y_max])
            # 更新心电数据
            if plot_first == True:
                ax.plot(np.arange(0,100)/fs, data, color='blue', linestyle='-', linewidth=1.5, label='ECG signal')
                plot_first = False
            else:
                ax.plot(np.arange(ii-100-n_overlap, ii)/fs,data[-100-n_overlap:], color='blue', linestyle='-', linewidth=1.5)

        if ii < 1600:
            text = 'Calibrating...'
        else:
            # 检测R波并标注
            if peak_first == True:
                # 由于第一个峰缺乏先验判断，无法确认是否为T波，因此舍弃
                QRS, T1, T2, SPK, NPK = detect_R(data,0,fs,peak_first=peak_first,lc=lc, hc=hc, order=order, wl=wl, param_N=param_N, param_S=param_S, param_T=param_T, iter_peak=iter_peak, iter_NS=iter_NS, iter_re_peak=iter_re_peak, iter_re_NS=iter_re_NS, thr_T=thr_T, thr_T_quotient=thr_T_quotient)
                QRS = QRS[1:]
                QRS_recent = QRS

                peak_first = False
                RR = np.diff(QRS_recent)
                RR_mean = np.mean(RR)
                RR_new = RR[(RR>0.92*RR_mean) & (RR<1.16*RR_mean)]
                RR_new_mean = np.mean(RR_new)
                idx = np.where(RR > 1.66*RR_new_mean)[0]
                if len(idx) > 0:
                    for i in range(len(idx)):
                        t1 = QRS_recent[idx[i]]
                        t2 = QRS_recent[idx[i]+1]
                        # 开始复检
                        QRS_new, T1, T2, SPK, NPK = detect_R(data[t1-n_overlap:t2+n_overlap],QRS_recent[idx[i]],fs,peak_first=peak_first,retest=True, T1=T1, T2=T2, SPK=SPK, NPK=NPK,data_size_previous=t1-n_overlap, lc=lc, hc=hc, order=order, wl=wl, param_N=param_N, param_S=param_S, param_T=param_T, iter_peak=iter_peak, iter_NS=iter_NS, iter_re_peak=iter_re_peak, iter_re_NS=iter_re_NS, thr_T=thr_T, thr_T_quotient=thr_T_quotient)
                        if len(QRS_new) != 0:
                            QRS_new = QRS_new + t1
                            if abs(t1 - QRS_new[0]) < 0.2*fs:
                                QRS_new = QRS_new[1:]
                            if len(QRS_new) != 0:
                                if abs(t2 - QRS_new[-1]) < 0.2*fs:
                                    QRS_new = QRS_new[:-1]
                            if len(QRS_new) != 0:
                                QRS_recent = np.insert(QRS_recent, idx[i]+1, QRS_new)

                
                QRS_old = QRS_recent
                # fine tune
                QRS_old = np.array([np.argmax(np.abs(data[max(0, int(i)-20):min(len(data), int(i)+20)])) + max(0, int(i)-20) for i in QRS_old])
                QRS_to_draw = np.sort(QRS_old)[:-1]
                x_coords = QRS_to_draw/fs            
                # y_coords = [data[np.argmax(np.abs(data[max(0, int(i)-20):min(len(data), int(i)+20)])) + max(0, int(i)-20)] for i in QRS_to_draw]
                y_coords = data[QRS_to_draw]

                width_old = calculate_width(QRS_old)

                if len(QRS_old) >= 3:
                    for QRS_diagnose in range(1,len(QRS_old)-1):
                        diagnose(QRS_old[QRS_diagnose-1:QRS_diagnose+2],width_old[QRS_diagnose-1:QRS_diagnose+2],fs,thr_width,thr_quotient)
                bpm, RR, status = calculate_heart_rate(QRS_recent, fs)

                if visible:
                    ax.scatter(x_coords, y_coords, s=50, marker='o', label='peaks', color='orange')
                    for i in premature_new:
                        ax.scatter(i/fs, y_max-0.1, s=50, marker='^', color='red')
                    premature_new = []

                    fill_start = QRS_recent[0]
                    fill_end = QRS_recent[-1]
                    fill_x = np.arange(fill_start,fill_end+1)/fs
                    fill = plt.fill_between(fill_x, y_min * np.ones_like(fill_x), y_max * np.ones_like(fill_x), color=colors[status],alpha=0.2)
                
            else:
                QRS, T1, T2, SPK, NPK = detect_R(data[-chunk_size-n_overlap:],QRS_old[-1], fs, peak_first=False, T1=T1, T2=T2, SPK=SPK, NPK=NPK,data_size_previous=ii-chunk_size-n_overlap, lc=lc, hc=hc, order=order, wl=wl, param_N=param_N, param_S=param_S, param_T=param_T, iter_peak=iter_peak, iter_NS=iter_NS, iter_re_peak=iter_re_peak, iter_re_NS=iter_re_NS, thr_T=thr_T, thr_T_quotient=thr_T_quotient)
                if len(QRS) == 0:
                    continue
                QRS = QRS + ii - chunk_size - n_overlap
                if abs(QRS_old[-1] - QRS[0]) < 0.2*fs:         # 检测重复
                    QRS = QRS[1:]
                    if len(QRS) == 0:
                        continue
                
                if len(QRS) + len(QRS_old) <= QRS_for_bpm:
                    QRS_recent = np.concatenate((QRS_old,QRS))
                else:
                    QRS_recent = np.concatenate((QRS_old[-(QRS_for_bpm - len(QRS)):],QRS))

                bpm, RR, status = calculate_heart_rate(QRS_recent, fs)
                
                RR_mean = np.mean(RR)
                RR_new = RR[(RR>0.92*RR_mean) & (RR<1.16*RR_mean)]
                while len(RR_new) == 0:
                    # 去掉RR中距离平均值最大的值
                    RR = np.delete(RR, np.argmax(np.abs(RR-RR_mean)))
                    RR_mean = np.mean(RR)
                    RR_new = RR[(RR>0.92*RR_mean) & (RR<1.16*RR_mean)]
                RR_new_mean = np.mean(RR_new)
                # idx = np.where(RR > 1.66*RR_new_mean)[0]
                idx = np.where(RR > 1.56*RR_new_mean)[0]
                if len(idx) > 0:
                    for i in range(len(idx)):
                        t1 = int(QRS_recent[idx[i]])
                        t2 = int(QRS_recent[idx[i]+1])
                        # 开始复检
                        QRS_new, T1, T2, SPK, NPK = detect_R(data[t1-n_overlap:t2+n_overlap],QRS_recent[idx[i]],fs,peak_first=peak_first,T1=T1,T2=T2,SPK=SPK,NPK=NPK,data_size_previous=ii-chunk_size-n_overlap,retest=True, lc=lc, hc=hc, order=order, wl=wl, param_N=param_N, param_S=param_S, param_T=param_T, iter_peak=iter_peak, iter_NS=iter_NS, iter_re_peak=iter_re_peak, iter_re_NS=iter_re_NS, thr_T=thr_T, thr_T_quotient=thr_T_quotient)
                        if len(QRS_new) != 0:
                            QRS_new = QRS_new + t1
                            if abs(t1 - QRS_new[0]) < 0.2*fs:         # 检测重复
                                QRS_new = QRS_new[1:]
                            if len(QRS_new) != 0:
                                if abs(t2 - QRS_new[-1]) < 0.2*fs:         # 检测重复
                                    QRS_new = QRS_new[:-1]
                            if len(QRS_new) != 0:
                                QRS_recent = np.insert(QRS_recent, idx[i]+1, QRS_new)
                                QRS_new_all = np.append(QRS_new_all,QRS_new)
                                bpm, RR, status = calculate_heart_rate(QRS_recent[len(QRS_new_all):], fs)
                            
                QRS_update = np.concatenate((QRS,QRS_new_all,[QRS_old[-1]]))
                QRS_update = np.array(QRS_update).astype(int)
                QRS_update = np.array([np.argmax(np.abs(data[max(0, int(i)-20):min(len(data), int(i)+20)])) + max(0, int(i)-20) for i in QRS_update])
                width_update = calculate_width(QRS_update)

                QRS_old = QRS_old[:-1]
                width_old = width_old[:-1]

                QRS_old = np.concatenate((QRS_update, QRS_old))
                width_old = np.concatenate((width_update, width_old))

                sort_indices = np.argsort(QRS_old)
                QRS_old = QRS_old[sort_indices]
                width_old = width_old[sort_indices]
                
                QRS_to_draw = np.sort(QRS_update)[:-1]
                x_coords = QRS_to_draw/fs
                y_coords = data[QRS_to_draw]

                indices = np.where(np.in1d(QRS_old, QRS_update))[0]
                for i in indices:
                    diagnose(QRS_old[i-9:i+2], width_old[i-9:i+2],fs,thr_width,thr_quotient)

                if visible:
                    ax.scatter(x_coords, y_coords, s=50, marker='o', color='orange')
                    for i in premature_new:
                        ax.scatter(i/fs, y_max-0.1, s=50, marker='^', color='red')
                    premature_new = []

                    fill_start = QRS_recent[0]
                    fill_end = QRS_recent[-1]
                    fill_x = np.arange(fill_start,fill_end+1)/fs
                    fill.remove()
                    fill = plt.fill_between(fill_x, y_min * np.ones_like(fill_x), y_max * np.ones_like(fill_x), color=colors[status],alpha=0.2)

            if visible:
                text = f'Heart rate: {bpm:.2f} bpm'
                ax.legend(loc=[0.9, 1.03])


        if visible:
            for txt in ax.texts:
                txt.remove()
            ax.text(0.93, 0.98, text, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

            plt.draw()
            plt.pause(0.00001)

    bpm, RR, status = calculate_heart_rate(QRS_old, fs)
    text = f"finish, {bpm:.2f} bpm"

    if visible:
        QRS_to_draw = np.argmax(np.abs(data[max(0, int(QRS_old[-1])-20):min(len(data), int(QRS_old[-1])+20)])) + max(0, int(QRS_old[-1])-20)
        x_coords = QRS_to_draw/fs
        y_coords = data[QRS_to_draw]
        ax.scatter(x_coords, y_coords, s=50, marker='o', color='orange')
        
        for txt in ax.texts:
            txt.remove()
        fill.remove()

        ax.text(0.95, 0.98, text, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    if visible or show_result:
        if Bradycardia:
            Bradycardia.merge_overlaps()
            for i in Bradycardia:
                Bradycardia_range = np.arange(i.begin, i.end+1)
                if Normal:
                    Normal.chop(i.begin, i.end)
                ax.fill_between(Bradycardia_range/fs, y_min * np.ones_like(Bradycardia_range), y_max * np.ones_like(Bradycardia_range), color=colors["Bradycardia"],alpha=0.2)
            

        if Tachycardia:
            Tachycardia.merge_overlaps()
            for i in Tachycardia:
                Tachycardia_range = np.arange(i.begin, i.end+1)
                if Normal:
                    Normal.chop(i.begin, i.end)
                ax.fill_between(Tachycardia_range/fs, y_min * np.ones_like(Tachycardia_range), y_max * np.ones_like(Tachycardia_range), color=colors["Tachycardia"],alpha=0.2)

        if Normal:
            Normal.merge_overlaps()
            for i in Normal:
                Normal_range = np.arange(i.begin, i.end+1)
                ax.fill_between(Normal_range/fs, y_min * np.ones_like(Normal_range), y_max * np.ones_like(Normal_range), color=colors["Normal"],alpha=0.2)

    predict = QRS_old / fs

    if evaluate:
        truth,V_truth,A_truth = read_annot(data_idx)
        premature_truth = np.concatenate((V_truth,A_truth))
        tp, fp, fn = compare_annot(truth, predict)
        tp_all[data_idx-1] = tp
        fp_all[data_idx-1] = fp
        fn_all[data_idx-1] = fn
        premature_predict = np.array(premature)/fs
        premature_tp, premature_fp, premature_fn = compare_annot(premature_truth, premature_predict)
        premature_tp_all[data_idx-1] = premature_tp
        premature_fp_all[data_idx-1] = premature_fp
        premature_fn_all[data_idx-1] = premature_fn

        accuracy[data_idx-1] = np.divide(np.float64(tp), np.float64(tp + fp), out=np.full_like(tp, np.nan, dtype=np.float64), where=(tp+fp)!=0)
        sensitivity[data_idx-1] = np.divide(np.float64(tp), np.float64(tp + fn), out=np.full_like(tp, np.nan, dtype=np.float64), where=(tp+fn)!=0)
        premature_accuracy[data_idx-1] = np.divide(np.float64(premature_tp), np.float64(premature_tp + premature_fp), out=np.full_like(premature_tp, np.nan, dtype=np.float64), where=(premature_tp+premature_fp)!=0)
        premature_sensitivity[data_idx-1] = np.divide(np.float64(premature_tp), np.float64(premature_tp + premature_fn), out=np.full_like(premature_tp, np.nan, dtype=np.float64), where=(premature_tp+premature_fn)!=0)

        f1[data_idx-1] = 2 * accuracy[data_idx-1] * sensitivity[data_idx-1] / (accuracy[data_idx-1] + sensitivity[data_idx-1])
        f1_premature[data_idx-1] = 2 * premature_accuracy[data_idx-1] * premature_sensitivity[data_idx-1] / (premature_accuracy[data_idx-1] + premature_sensitivity[data_idx-1])


    if show_result:
        ax.set_xlim([0, len(data_all)/fs])
        ax.set_ylim([y_min, y_max])
        ax.plot(np.arange(0,len(data_all))/fs, data_all, color='blue', linestyle='-', linewidth=1.5, label='ECG signal')
        ax.scatter(QRS_old/fs, data_all[QRS_old], s=50, marker='o', color='orange')
        for i in premature:
            ax.scatter(i/fs, y_max-0.1, s=50, marker='^', color='red')
        ax.text(0.95, 0.98, text, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        plt.savefig(os.path.join(result_path, f'{data_idx}.png'))
        plt.show(block=False)
        plt.pause(5)
        plt.close()

        

    if visible:
        plt.show(block=False)
        plt.pause(1)
        plt.close()

if evaluate:
    tp_all[10] = np.sum(tp_all)
    fp_all[10] = np.sum(fp_all)
    fn_all[10] = np.sum(fn_all)
    tp_all[11] = tp_all[10] - tp_all[3]
    fp_all[11] = fp_all[10] - fp_all[3]
    fn_all[11] = fn_all[10] - fn_all[3]

    premature_tp_all[10] = np.sum(premature_tp_all)
    premature_fp_all[10] = np.sum(premature_fp_all)
    premature_fn_all[10] = np.sum(premature_fn_all)
    premature_tp_all[11] = premature_tp_all[10] - premature_tp_all[3]
    premature_fp_all[11] = premature_fp_all[10] - premature_fp_all[3]
    premature_fn_all[11] = premature_fn_all[10] - premature_fn_all[3]

    accuracy = np.divide(tp_all, tp_all + fp_all, out=np.full_like(tp_all, np.nan, dtype=np.float64), where=(tp_all+fp_all)!=0)
    sensitivity = np.divide(tp_all, tp_all + fn_all, out=np.full_like(tp_all, np.nan, dtype=np.float64), where=(tp_all+fn_all)!=0)
    premature_accuracy = np.divide(premature_tp_all, premature_tp_all + premature_fp_all, out=np.full_like(premature_tp_all, np.nan, dtype=np.float64), where=(premature_tp_all+premature_fp_all)!=0)
    premature_sensitivity = np.divide(premature_tp_all, premature_tp_all + premature_fn_all, out=np.full_like(premature_tp_all, np.nan, dtype=np.float64), where=(premature_tp_all+premature_fn_all)!=0)

    f1 = 2 * accuracy * sensitivity / (accuracy + sensitivity)
    f1_premature = 2 * premature_accuracy * premature_sensitivity / (premature_accuracy + premature_sensitivity)

    # 存入excel
    df = pd.DataFrame({
        'idx':np.concatenate((np.arange(1,11).astype(str),['all','outlier_free'])),
        'TP': tp_all,
        'FP': fp_all,
        'FN': fn_all,
        'Accuracy': accuracy,
        'Sensitivity': sensitivity,
        'F1': f1,
        'Premature TP': premature_tp_all,
        'Premature FP': premature_fp_all,
        'Premature FN': premature_fn_all,
        'Premature Accuracy': premature_accuracy,
        'Premature Sensitivity': premature_sensitivity,
        'Premature F1': f1_premature
    })

    df.to_excel(eval_result_path, index=False)