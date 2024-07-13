import numpy as np
from scipy.signal import butter, filtfilt, find_peaks

def detect_R(data_original,QRS_previous,fs=200,peak_first=True,T1=0,T2=0,SPK=0,NPK=0,data_size_previous=1600,retest=False,lc=5,hc=15,order=6,wl=30,param_N=0.75,param_S=0.25,param_T=0.5,iter_peak=0.125,iter_NS=0.875,iter_re_peak=0.25,iter_re_NS=0.75,thr_T=75,thr_T_quotient=0.25):
    # 高通滤波
    lowcut = lc/fs
    b, a = butter(order, lowcut , btype='high')
    high_passed = filtfilt(b, a, data_original)

    # 低通滤波
    highcut = hc/fs
    b, a = butter(order, highcut , btype='low')
    low_passed = filtfilt(b, a, high_passed)

    # 差分并平方
    diff = np.diff(low_passed)
    diff = np.insert(diff, 0, 0)
    diff = diff * fs
    diff_square = diff * diff

    # 加窗平均
    window = np.ones(wl) / wl
    data_average = np.convolve(diff_square, window, 'same')

    # 寻找峰值
    peaks, _ = find_peaks(data_average, prominence=1, distance=1)

    # 阈值迭代
    QRS = []

    # 赋初值
    if peak_first:      
        reshaped_data = data_average[:1600].reshape(-1, 200)
        max_values = np.max(reshaped_data, axis=1)
        SPK = np.mean(max_values)
        NPK = 0
        T1 = param_N*NPK + param_S*SPK
        T2 = param_T*T1

    pass_label = True
    for i in range(len(peaks)):
        if len(QRS) != 0:
            QRS_last = QRS[-1]
        else:
            QRS_last = QRS_previous - data_size_previous

        if pass_label:
            if (peaks[i] - QRS_last < 0.2*fs):
                continue
            else:
                pass_label = False

        threshold = T1*(1-retest) + T2*retest

        if data_average[peaks[i]] >= threshold:
            pass_label = True
            SPK_previous = SPK
            if retest == 0:
                SPK = iter_peak * data_average[peaks[i]] + iter_NS * SPK_previous
            else:
                SPK = iter_re_peak * data_average[peaks[i]] + iter_re_NS * SPK_previous
            QRS = np.append(QRS, peaks[i])

        else:
            NPK_previous = NPK
            NPK = iter_peak * data_average[peaks[i]] + iter_NS * NPK_previous
        
        T1 = param_N*NPK + param_S*SPK
        T2 = param_T*T1

    QRS = np.array(QRS).astype(int)


    # 剔除T波
    RR = np.diff(np.concatenate(([QRS_previous-data_size_previous],QRS)))
    idx = np.where(RR<thr_T)[0]
    T_suspect = QRS[idx]
    slope_T_suspect = data_average[T_suspect]
    slope_previous_R = data_average[T_suspect - 1] 
    quotient =  slope_previous_R / slope_T_suspect
    T = T_suspect[quotient > thr_T_quotient]
    QRS_new = np.setdiff1d(QRS, T)


    return QRS_new, T1, T2, SPK, NPK




