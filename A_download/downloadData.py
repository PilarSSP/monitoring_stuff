"""
downloadData.py

This script automates the download of continuous seismic waveform data from FDSN-compatible web services
(e.g., IGN, IRIS), using a list of stations and optional date ranges as input. It includes parallel downloading,
resampling, response removal, quality control, and SAC output formatting.

Usage:
    - Download data for today:
        python downloadData.py
    - Download data for a specific day:
        python downloadData.py 2023-01-31
    - Download data for a date range:
        python downloadData.py 2023-01-01 2023-01-31

Output:
    - Waveforms are saved in SAC format: `NET.STA.CHN.YYYY.JJJ.HH.MM.SS.SACv`
    - Stored under: `output_dir/NET/STA/CHN/YYYY/JJJ/`
    - Log files include: Maximum amplitude and standard deviation of the day-long and dh-long traces

Configuration:
    All input/output paths, filter settings, and client options are specified in config.py.

Author:
    Pilar Sánchez-Pastor, 2025 (psanchez@geo3bcn.csic.es)
"""

# -----------------   HEADER   ------------------
import os
import obspy as obs
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
import time
import numpy as np
import pandas as pd
from pathlib import Path
import sys

# Determine the project root dynamically (two levels up from this script)
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from config import CONFIG  # Import configuration dictionary

# Processing and filter settings
bp = CONFIG['bandpass_filter']                      # Bandpass filter for instrument response deconvolution
freq_new = CONFIG['resample_rate']                  # Target sampling rate after resampling (in Hz)

# Directory structure
dh = CONFIG['window_hours']                          # Length of each data segment in hours (e.g., 0.25 = 15 min)
dh_flag = CONFIG['window_label']                     # Label used for organizing output folders by window length
dir_out0 = CONFIG['output_base_dir']                 # Base output directory (from config)
dir_out = f"{dir_out0}/SACv{dh_flag}"                    # Final output path for SAC files

# Input file paths
in_file = project_root / "files" / CONFIG['station_list_file']    # File with list of stations to download
st_info = project_root / "files" / CONFIG['station_info_file']    # File with metadata for all stations

# Channel selection
ch = CONFIG['component']                             # Component to download (e.g., Z for vertical)
channel_priority = CONFIG['channel_priority']        # Channel priority order if multiple sensors are available for the same station

# Data source
fdsn_client = Client(CONFIG['fdsn_url'])             # FDSN web service client for waveform data download

max_workers = CONFIG.get('max_workers', 3)  # fallback = 3

# --------------------- Download Function ----------------------

def download_data(date_str: str, network, station, channel="*H"+ch):
    """
    Download seismic data for a station and save as SAC.
    Errors are caught to avoid breaking the process.
    """
    try:
        # Extract station coordinates
        with open(st_info, 'r') as info0:
            lines = info0.readlines()
        in_net, in_stat, in_lat, in_lon = zip(*[x.split()[:4] for x in lines])

        # Create a text file for trace information
        trace_dir = Path(dir_out0) / 'traces'
        trace_dir.mkdir(parents=True, exist_ok=True)
        tr_file = trace_dir / f"{network}_{station}_{ch}_{dh_flag}.txt"
        is_new_file = not tr_file.exists()

        with open(tr_file, "a") as fi:
            if is_new_file:
                fi.write("Trace\tAmax\tstd_day\tstd_tr\n")

            # Define time range
            date = UTCDateTime(date_str)
            t1 = UTCDateTime(date.year, date.month, date.day, 0, 0, 0)
            t2 = t1 + 24*3600

            # Downloading waveform and instrument response
            st = fdsn_client.get_waveforms(network=network, station=station, location="*", channel=channel, starttime=t1, endtime=t2, attach_response=True)
            
            # Check if multiple channels exist and keep only the best
            if len(set(tr.stats.channel for tr in st)) > 1:
                for ch0 in channel_priority:
                    selected = st.select(channel=ch0)
                    if len(selected) > 0:
                        st = selected
                        print(f"[INFO] Selected channel {ch0} for {st[0].id}")
                        break

            st.detrend(type='constant')  # rmean
            st.detrend(type='linear')  # rtrend
            freq = st[0].stats.sampling_rate
            st.merge(method=1, fill_value=0)  # merge traces

            year = f"{t1.year:04d}"
            day = f"{t1.julday:03d}"
            dir_outtr = f"{dir_out}/{network}/{station}/{ch}/{year}/{day}"
            os.makedirs(dir_outtr, exist_ok=True)

            # Consider traces longer than 1 hour
            if st[0].stats.npts < 3600 * freq:
                print(f"[SKIPPED] {network}.{station} on {date_str} → Short trace")
                return
                
            # Skip trace if its native sampling rate is lower than the desired one
            if freq < freq_new:
                print(f"[SKIPPED] {network}.{station} on {date_str}: Sampling freq {freq} Hz")
                return

            # Trim to full hours or minutes
            if st[0].stats.starttime <= t1:
                st.trim(t1, t2, nearest_sample=False)
            else:
                t1_0 = st[0].stats.starttime
                r = int(t1_0.hour % dh)
                if (dh - r) == 0:
                    t1_1 = obs.core.utcdatetime.UTCDateTime(year=t1.year, julday=t1.julday, hour=t1_0.hour, minute=1)
                else:
                    t1_1 = obs.core.utcdatetime.UTCDateTime(year=t1.year, julday=t1.julday, hour=t1_0.hour+int((dh - r)))
                st.trim(t1_1, t2, nearest_sample=False)
                del (t1_0, t1_1)
                
            # Downsampling if needed
            if freq_new != freq:
                st.filter('lowpass', freq=freq_new / 2, zerophase=True)
                if st[0].stats.starttime <= t1:
                    st.interpolate(sampling_rate=freq_new, starttime=t1)
                    
            # Removing instrument response
            st.remove_response(output='VEL', pre_filt=bp, zero_mean=True, taper=False)  # Deconvolve response

            # Add station location to SAC header
            tr = st[0]
            tr.stats.sac = obs.core.AttribDict()
            if station not in in_stat:
                print(f"[WARN] Station {station} not found in metadata.")
                return
            pos = in_stat.index(station)
            tr.stats.sac.stla = float(in_lat[pos])
            tr.stats.sac.stlo = float(in_lon[pos])

            # Standard deviation over the entire day-long trace
            std_day = tr.std()

            # Trim in dh hours
            tc1 = t1
            tc2 = t1 + 3600 * dh - 1e-6
            to = tr.stats.starttime
            te = tr.stats.endtime

            while tc2 <= te + 60:
                if tc2 <= to:
                    tc1 = tc2 + 1e-6
                    tc2 = tc1 + 3600 * dh - 1e-6
                    continue

                tr1 = tr.copy()
                tr1.trim(tc1, tc2, nearest_sample=False)
                tt = tr1.stats.starttime
                ti = "%.2d" % tt.hour + '.' + "%.2d" % tt.minute + '.' + "%.2d" % tt.second
                trace = network + '.' + station + '.' + tr1.stats.channel + '.' + year + '.' + day + '.' + ti + '.SACv'
                file_out = dir_outtr + '/' + trace
                tr1.write(file_out, format='SAC')

                # Maximum amplitude and standard deviation of the dh-long trace
                Am = max(max(tr1), abs(min(tr1)))
                std_tr = tr1.std()

                fi.write("%s\t%.4e\t%.3e\t%.3e\n" % (trace, Am, std_day, std_tr))

                tc1 = tc2 + 1e-6
                tc2 = tc1 + 3600 * dh - 1e-6

                del (tr1, tt, ti, file_out)

            print(f"[OK]  {network}.{station}")

    except Exception as e:
        print(f"[ERROR] {network}.{station} on {date_str}: {e}")
        return

    # Sort and clean traces file
    data = pd.read_csv(tr_file, sep="\t", skiprows=[0], names=['Trace', 'Amax', 'std_day', 'std_tr'])
    data.sort_values(data.columns[0], inplace=True)
    data_ok = data.drop_duplicates()
    np.savetxt(tr_file, data_ok.values, fmt=('%s', '%.3e', '%.3e', '%.3e'), delimiter="\t",
               header=str(data.columns.values))


# --------------------- Main Execution ----------------------

if __name__ == "__main__":
    from concurrent.futures import ThreadPoolExecutor, as_completed
    from datetime import datetime, timedelta
    import sys
    
    time0 = time.time()

    # Date range input from command line or default to today
    date_start = sys.argv[1] if len(sys.argv) > 1 else UTCDateTime().date.isoformat()
    date_end = sys.argv[2] if len(sys.argv) > 2 else date_start

    start_date = datetime.fromisoformat(date_start)
    end_date = datetime.fromisoformat(date_end)

    # Read station list from input file (path now from config)
    if os.path.exists(in_file):
        df_in = pd.read_csv(in_file, sep='\s+', header=None, names=["network", "station"])
    else:
        raise FileNotFoundError(f"Input file {in_file} not found.")

    current_date = start_date
    while current_date <= end_date:
        date_str = current_date.strftime("%Y-%m-%d")
        print(f"\n=== Processing day: {date_str} ===")

        def worker(row):
            network = row["network"]
            station = row["station"]
            try:
                download_data(date_str, network=network, station=station, channel='*H'+ch)
            except Exception as e:
                print(f"[ERROR] processing {network}.{station} en {date_str}: {e}")

        # Run parallel downloads
        n_workers = min(max_workers, len(df_in))
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            futures = [executor.submit(worker, row) for idx, row in df_in.iterrows()]
            for future in as_completed(futures):
                future.result()

        current_date += timedelta(days=1)

    time1 = time.time()
    print('\n\n   ->   Completed process in ', time.strftime('%H:%M:%S', time.gmtime(time.time() - time0)))
    print("\n\n\n( )( )\n(^ .^)  < Done!\n(”__”)\n")
