import os
import os.path as osp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from pygam import s, f, ExpectileGAM
from functools import partial

from joblib import Parallel, delayed

import matplotlib.pyplot as plt
from tqdm import tqdm
import argparse

# Input file names for different strain types
INPUT_FNAMES = {
    "radial": "strain_sa_radial.csv",
    "circum": "strain_sa_circum.csv",
    "longit": "strain_la_4ch_longit.csv",
}

FAILED_DETECT = "failed.txt"
MISSING_FILE = "missing_input.txt"

def find_peak_gam(x, y, expectile=0.5, max_skip_right=1 / 4, invert_sign=False):
    """
    Find the first derivative peak from the input signal.
    The `max_skip_right` argument determines the right proportion of the signal that is
    initially discarded to avoid the last peak.
    """
    # --------------------------------------------- Find the signal peak --------------------------------------------- #
    max_skip = int(len(y) * max_skip_right)
    crop_left = 5
    crop_right = 5

    # Sometimes there's an artefact peak at the beginning of the signal. We cut it here
    x_crop = x[crop_left:-crop_right]

    gam = ExpectileGAM(expectile=expectile)
    gam = gam.gridsearch(
        x.reshape(-1, 1),
        y.reshape(
            -1,
        ),
        progress=False,
    )
    y_gam_crop = gam.predict(x_crop.reshape(-1, 1))

    # Determine if the signal must be inverted. This is done by looking at the first high peak which should correspond
    # to the main peak of the signal. If this is negative, then the signal must be inverted before processing.
    gam_peaks_pos, gam_props_pos = find_peaks(y_gam_crop, rel_height=0.5, width=0)
    gam_peaks_neg, gam_props_neg = find_peaks(-y_gam_crop, rel_height=0.5, width=0)
    if len(gam_peaks_pos) == 0 and len(gam_peaks_neg) == 0:
        raise RuntimeError("No signal peak found.")
    elif len(gam_peaks_pos) == 0:
        invert_sign = True
    elif len(gam_peaks_neg) == 0:
        invert_sign = False
    elif np.min(gam_peaks_neg) < np.min(gam_peaks_pos):
        invert_sign = True
    else:
        invert_sign = False

    # Smooth the whole signal
    y_gam = gam.predict(x.reshape(-1, 1))

    # Here we invert the signal, if necessary
    if invert_sign:
        y = -y
        y_gam = -y_gam
        gam_peaks = gam_peaks_neg
        gam_props = gam_props_neg
    else:
        gam_peaks = gam_peaks_pos
        gam_props = gam_props_pos

    # As we are processing the whole signal now, we need to shift the peaks found in the cropped signal
    gam_peaks += crop_left
    gam_props["left_ips"] += crop_left
    gam_props["right_ips"] += crop_left

    der_ = np.diff(y)
    der_gam = np.diff(y_gam)

    # Find closest real peak to the one from the smoothed signal
    opt_ind = np.argmax(y[gam_peaks])
    left_ip = int(np.floor(gam_props["left_ips"][opt_ind]))
    right_ip = int(np.ceil(gam_props["right_ips"][opt_ind]))

    # The candidate peaks are searched in the left_ips, right_ips interval.
    s_peaks = find_peaks(y)[0]
    cand = np.asarray([int(x) for x in s_peaks if x >= left_ip and x <= right_ip])
    # If no candidates are found in the interval, take the closest peak
    if len(cand) == 0:
        S_start = s_peaks[np.argmin(np.abs(s_peaks - gam_peaks[opt_ind]))]
    else:
        S_start = cand[np.argmax(y[cand])]

    # Now search the derivative peak, iteratively. Start by applying the max crop to the right and then reduce it
    # until it finds a high negative peak in the derivative.
    found = False
    for sp in range(max_skip, -1, -1):
        try:
            query = der_gam[S_start:-sp]
            inds, props = find_peaks(-query, rel_height=0.75, width=0)
            if len(inds) == 0:
                raise RuntimeError(f"skip {sp}: no der. peak")
            else:
                opt_ind = np.argmin(query[inds])
                peak_ind = inds[opt_ind] + S_start
                left_ip = int(np.floor(props["left_ips"][opt_ind])) + S_start
                right_ip = int(np.ceil(props["right_ips"][opt_ind])) + S_start
                found = True
        except:
            continue
        if found:
            break
    if not found:
        raise RuntimeError("No der. peak found")

    # Determine the closest lowest peak in the derivative
    der_peaks = find_peaks(-der_)[0]
    cand = np.asarray([int(x) for x in der_peaks if left_ip <= x <= right_ip])
    if len(cand) != 0:
        der_peak_ = cand[np.argmin(der_[cand])]
    else:
        # If no peaks are found in the interval, take the closest.
        der_peak_ = der_peaks[np.argmin(np.abs(der_peaks - peak_ind))]

    if invert_sign:
        der_ = -der_

    return S_start, der_peak_, der_


def plot_results(tpoints, signal, der, sig_peak, der_peak, signal_name, subject, plots_dir):
    os.makedirs(os.path.join(plots_dir, subject), exist_ok=True)

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(tpoints, signal, label="signal", c="red")
    ax2.plot(tpoints[:-1], der, label="der", c="green")
    if sig_peak is not None:
        ax1.plot(tpoints[sig_peak], signal[sig_peak], "ko", label="signal peak")
    if der_peak is not None:
        ax2.plot(tpoints[der_peak], der[der_peak], "kv", label="der peak")

    ax1.set_ylabel("Signal")
    ax2.set_ylabel("Derivative")

    ax2.legend(loc="upper right")
    ax1.legend(loc="upper left")
    ax1.set_title(f"Signal: {signal_name}")
    fig.savefig(osp.join(plots_dir, subject, f"{subject}_{signal_name}.png"))
    plt.close()

def merge_results(input_dir, save_path):
    list_csv = os.listdir(input_dir)
    subjects = [x.split("_")[0] for x in list_csv]
    merged_df = []
    for f in tqdm(list_csv):
        values = pd.read_csv(osp.join(input_dir, f), index_col=0)
        merged_df.append(values.to_numpy().squeeze())
    merged_df = np.r_[merged_df]
    merged_df = pd.DataFrame(data=merged_df, columns=values.columns)
    merged_df["ID"] = subjects
    merged_df.to_csv(save_path)
    return None

def parse_args():
    parser = argparse.ArgumentParser(description="Find peaks strain UKBB")
    parser.add_argument("signal", type=str, choices=["radial", "circum", "longit"], default="radial")
    parser.add_argument("input_dir", type=str)
    parser.add_argument("output_dir", type=str)
    parser.add_argument("--only_global", action="store_true", dest="only_global", default=False)
    parser.add_argument(
        "--list_subjects",
        dest="list_subjects",
        type=str,
        help="Text file containing the list of IDs to analyse (txt format).",
        default=None,
    )
    parser.add_argument("--do_plots", dest="do_plots", action="store_true", default=False)

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    def _process_subjects(subject, only_global: bool, invert_signal: bool):
        if osp.isfile(osp.join(args.input_dir, subject)):
            return None

        subject_dir = osp.join(args.input_dir, subject)

        # if osp.isfile(osp.join(save_dir, f"{subject}.csv")):
        #     return None
        try:
            visit_folders = os.listdir(subject_dir)
        except Exception as e:
            print(f"Error accessing subject directory {subject_dir}: {e}")
            with open(osp.join(out_dir, MISSING_FILE), "a") as f:
                f.writelines(f"{subject}\n")
            return None     #subject_dir  = subject_dir + "/visit_2/"
        for visit_folder in visit_folders:
            try:
                visit_dir = osp.join(subject_dir, visit_folder)
                input_file_path = osp.join(visit_dir, input_fname)

                if not osp.isfile(input_file_path):
                    with open(osp.join(out_dir, MISSING_FILE), "a") as f:
                        f.writelines(f"{subject}\n")
                    continue

                signals = pd.read_csv(osp.join(subject_dir, input_fname), index_col=0)
                # Remove all-zero columns
                signals = signals.loc[:, (signals != 0).any(axis=0)]

                # Time points from column names
                tpoints = np.asarray(signals.columns, dtype=float)
                tpoints = np.round(tpoints, 4)  # Round the time points

                # Take signal names from rownames

                # Prepare output
                results_dict = {
                    "signal_name": [],
                    "sig_peak_idx": [],
                    "sig_peak_time": [],
                    "sig_peak_value": [],
                    "der_peak_idx": [],
                    "der_peak_time": [],
                    "der_peak_value": [],
                }

                signal_names = np.asarray([str(x) for x in signals.index])

                if only_global:
                    # Use only Global
                    signal_names = ["Global"]

                for i in range(len(signal_names)):
                    signal = signals.loc[signals.index == signal_names[i], :].to_numpy().squeeze()

                    try:
                        sig_peak, der_peak, der = find_peak_gam(tpoints, signal.copy(), invert_sign=invert_signal)
                    except:
                        with open(osp.join(out_dir, FAILED_DETECT), "a") as f:
                            f.writelines(f"{subject}, {signal_names[i]}\n")
                        der = np.diff(signal)
                        sig_peak = None
                        der_peak = None
                        results_dict["signal_name"].append(signal_names[i])
                        results_dict["sig_peak_idx"].append("-")
                        results_dict["sig_peak_time"].append("-")
                        results_dict["sig_peak_value"].append("-")
                        results_dict["der_peak_idx"].append("-")
                        results_dict["der_peak_time"].append("-")
                        results_dict["der_peak_value"].append("-")
                        plot_results(tpoints, signal, der, sig_peak, der_peak, signal_names[i], subject, plots_failed_dir)
                        continue

                    if args.do_plots:
                        plot_results(tpoints, signal, der, sig_peak, der_peak, signal_names[i], subject, plots_dir)

                    results_dict["signal_name"].append(signal_names[i])
                    results_dict["sig_peak_idx"].append(sig_peak)
                    results_dict["sig_peak_time"].append(tpoints[sig_peak])
                    results_dict["sig_peak_value"].append(signal[sig_peak])
                    results_dict["der_peak_idx"].append(der_peak)
                    results_dict["der_peak_time"].append(tpoints[der_peak])
                    results_dict["der_peak_value"].append(der[der_peak])

                df = pd.DataFrame.from_dict(results_dict)
                df.to_csv(osp.join(save_dir, f"{subject}.csv"))
            except Exception as e:
                print(f"Error processing visit folder {visit_folder} for subject {subject}: {e}")
                continue


    args = parse_args()

    if args.list_subjects is None:
        subjects = [x for x in os.listdir(args.input_dir)]
        print(f"Found {len(subjects)} directories")
    else:
        subjects = np.genfromtxt(args.list_subjects, dtype=str)
        subjects = np.asarray(subjects).reshape(
            -1,
        )
        print(subjects)
        print(len(subjects))
        print(f"Processing {len(subjects)} subjects...")

    os.makedirs(args.output_dir, exist_ok=True)
    out_dir = osp.join(args.output_dir, args.signal)
    save_dir = osp.join(out_dir, "results")
    plots_dir = osp.join(out_dir, "plots")
    plots_failed_dir = osp.join(out_dir, "plots_failed")
    os.makedirs(save_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)
    os.makedirs(plots_failed_dir, exist_ok=True)

    if args.signal in ["longit", "circum"]:
        invert_signal = True
    else:
        invert_signal = False

    input_fname = INPUT_FNAMES[args.signal]

    Parallel(n_jobs=-1)(
        delayed(partial(_process_subjects, only_global=args.only_global, invert_signal=invert_signal))(subject)
        for subject in tqdm(sorted(subjects))
    )
    for subject in tqdm(subjects):
        _process_subjects(only_global=args.only_global, invert_signal=invert_signal, subject=subject)

    for subject in tqdm(subjects):
    
        if osp.isfile(osp.join(args.input_dir, subject)):
            continue
    
        subject_dir = osp.join(args.input_dir, subject)
    
        if osp.isfile(osp.join(save_dir, f"{subject}.csv")):
            continue
    
        if not osp.isfile(osp.join(subject_dir, input_fname)):
            with open(osp.join(out_dir, MISSING_FILE), 'a') as f:
                f.writelines(f'{subject}\n')
            continue
    
        signals = pd.read_csv(osp.join(subject_dir, input_fname),
            index_col=0)
        # Remove all-zero columns
        signals = signals.loc[:, (signals != 0).any(axis=0)]
    
        # Time points from column names
        tpoints = np.asarray(signals.columns, dtype=float)
        tpoints = np.round(tpoints, 4)  # Round the time points
    
        # Take signal names from rownames
    
        # # Use only Global
        # i = int(np.where(signal_names == 'Global')[0])
    
        # Prepare output
        results_dict = {
            'signal_name': [],
            'sig_peak_idx': [],
            'sig_peak_time': [],
            'sig_peak_value': [],
            'der_peak_idx': [],
            'der_peak_time': [],
            'der_peak_value': []
        }
    
        signal_names = np.asarray([str(x) for x in signals.index])
    
        for i in range(len(signal_names)):
            signal = signals.iloc[i, :].to_numpy().squeeze()
    
            try:
                sig_peak, der_peak, der = find_peak_gam(tpoints, signal.copy(), invert_sign=invert_signal)
                if invert_signal:
                    der = -der
            except:
                with open(osp.join(out_dir, FAILED_DETECT), 'a') as f:
                    f.writelines(f'{subject}, {signal_names[i]}\n')
                der = np.diff(signal)
                sig_peak = None
                der_peak = None
                results_dict['signal_name'].append(signal_names[i])
                results_dict['sig_peak_idx'].append("-")
                results_dict['sig_peak_time'].append("-")
                results_dict['sig_peak_value'].append("-")
                results_dict['der_peak_idx'].append("-")
                results_dict['der_peak_time'].append("-")
                results_dict['der_peak_value'].append("-")
                plot_results(tpoints, signal, der, sig_peak, der_peak, signal_names[i], subject, plots_failed_dir)
                continue
    
            plot_results(tpoints, signal, der, sig_peak, der_peak, signal_names[i], subject, plots_dir)
    
            results_dict['signal_name'].append(signal_names[i])
            results_dict['sig_peak_idx'].append(sig_peak)
            results_dict['sig_peak_time'].append(tpoints[sig_peak])
            results_dict['sig_peak_value'].append(signal[sig_peak])
            results_dict['der_peak_idx'].append(der_peak)
            results_dict['der_peak_time'].append(tpoints[der_peak])
            results_dict['der_peak_value'].append(der[der_peak])
    
        df = pd.DataFrame.from_dict(results_dict)
        df.to_csv(osp.join(save_dir, f'{subject}.csv'))

