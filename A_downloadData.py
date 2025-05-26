# This script downloads seismic data from the EIDA FDSN web service

# Created by Pilar Sánchez, 2025

# -----------------   HEADER   ------------------
import os
import obspy as obs
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
import time
import numpy as np
import pandas as pd
from pathlib import Path

bp = (0.01, 0.05, 10, 20)
freq_new = 20 # New sampling rate

dir0 = '/home/psanchez/SHmax/SIRENS'
dh = 1 # Units: hours. Example: dh=0.25 means traces of 15-minutes duration
dh_flag = '1h'
dir_out0 = '/dataorion/psanchez/sirens/data'
dir_out =  dir_out0 + '/SACv' + dh_flag # Output directory

in_file = dir0 + '/files/in_stations.txt'
st_info = dir0 + '/files/stations_info.txt'
ch = 'Z'
channel_priority = ["HHZ", "BHZ", "LHZ"]

fdsn_client = Client("http://service.iris.edu/")

time0 = time.time()

# --------------------- Download Function ----------------------

def download_data(date_str: str, network, station, channel="*H"+ch):
    """
    Download seismic data for a station and save as SAC.
    Errors are caught to avoid breaking the process.
    """
    try:
        # Read station info to get coordinates
        with open(st_info, 'r') as info0:
            lines = info0.readlines()[1:]
        in_net, in_stat, in_lat, in_lon = zip(*[x.split()[:4] for x in lines])

        # Prepare trace file for writing trace info
        trace_dir = Path(dir_out0) / 'traces'
        trace_dir.mkdir(parents=True, exist_ok=True)
        tr_file = trace_dir / f"{network}_{station}_{ch}_{dh_flag}.txt"
        is_new_file = not tr_file.exists()

        with open(tr_file, "a") as fi:
            if is_new_file:
                fi.write("Trace\tSamp_gap\tAmax\tstd_day\tstd_tr\n")

            # Define time range
            date = UTCDateTime(date_str)
            t1 = UTCDateTime(date.year, date.month, date.day, 0, 0, 0)
            t2 = t1 + 24*3600

            # Try downloading waveform
            st = fdsn_client.get_waveforms(network=network, station=station, location="*", channel=channel, starttime=t1, endtime=t2, attach_response=True)
            # Check if multiple vertical channels exist
            if len(set(tr.stats.channel for tr in st)) > 1:
                for ch0 in channel_priority:
                    selected = st.select(channel=ch0)
                    if len(selected) > 0:
                        st = selected  # keep only the best available channel
                        print(f"[INFO] Selected channel {ch0} for {st[0].id}")
                        break
            else:
                print(f"[INFO] Only one vertical channel found: {st[0].stats.channel}")

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

            Fnyq = freq / 2
            if freq_new / 2 > Fnyq:
                print(f"[SKIPPED] {network}.{station} on {date_str}: Filter freq > Nyquist ({freq_new / 2}Hz >= {Fnyq}Hz)")
                return

            # Trim to full hours (00:00:00.000)
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

            if freq_new != freq:
                st.filter('lowpass', freq=freq_new / 2, zerophase=True)
                if st[0].stats.starttime <= t1:
                    st.interpolate(sampling_rate=freq_new, starttime=t1)

            st.remove_response(output='VEL', pre_filt=bp, zero_mean=True, taper=False)  # Deconvolve response

            # Add station location to SAC header
            tr = st[0]
            tr.stats.sac = obs.core.AttribDict()
            pos = in_stat.index(station)
            tr.stats.sac.stla = float(in_lat[pos])
            tr.stats.sac.stlo = float(in_lon[pos])

            # Standard deviation
            std_day = tr.std()

            # Cut in dh hours
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

                # Check for gaps: Adapt the function to the dataset!
                diff = np.abs(np.diff(tr1))
                thr = 1e-4 * np.std(diff)
                gaps0 = np.where(np.abs(np.diff(tr1)) < thr)[0]
                if len(gaps0) > 0:
                    gaps = len(gaps0)
                else:
                    gaps = 0

                # Creating txt files with gap and amplitude information
                # Maximum amplitude (absolute value)
                Am = max(max(tr1), abs(min(tr1)))

                # Standard deviation
                std_tr = tr1.std()

                fi.write("%s\t%d\t%.4e\t%.3e\t%.3e\n" % (trace, gaps, Am, std_day, std_tr))


                tc1 = tc2 + 1e-6
                tc2 = tc1 + 3600 * dh - 1e-6

                del (tr1, tt, ti, file_out)

            print(f"[OK]  {network}.{station}")

    except Exception as e:
        print(f"[ERROR] {network}.{station} on {date_str}: {e}")
        return

    # Sort and clean traces file
    data = pd.read_csv(tr_file, sep="\t", skiprows=[0], names=['Trace', 'Samp_gap', 'Amax', 'std_day', 'std_tr'])
    data.sort_values(data.columns[0], inplace=True)
    data_ok = data.drop_duplicates()
    np.savetxt(tr_file, data_ok.values, fmt=('%s', '%d', '%.3e', '%.3e', '%.3e'), delimiter="\t",
               header=str(data.columns.values))  # Overwrite


# --------------------- Main Execution ----------------------

if __name__ == "__main__":
    from concurrent.futures import ThreadPoolExecutor, as_completed
    from datetime import datetime, timedelta
    import sys

    date_start = sys.argv[1] if len(sys.argv) > 1 else UTCDateTime().date.isoformat()
    date_end = sys.argv[2] if len(sys.argv) > 2 else date_start

    start_date = datetime.fromisoformat(date_start)
    end_date = datetime.fromisoformat(date_end)

    # Read station list
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
                #print(f"Downloading data for {network}.{station}")
                download_data(date_str, network=network, station=station, channel='*H'+ch)
            except Exception as e:
                print(f"[ERROR] processing {network}.{station} en {date_str}: {e}")

        # Run parallel downloads
        n_workers = min(8, len(df_in))  # Hasta 8 estaciones en paralelo
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            futures = [executor.submit(worker, row) for idx, row in df_in.iterrows()]
            for future in as_completed(futures):
                future.result()

        current_date += timedelta(days=1)

    time1 = time.time()
    print('   ->   Completed process in ', time.strftime('%H:%M:%S', time.gmtime(time.time() - time0)))
    print("\n\n\n( )( )\n(^ .^)  < Done!\n(”__”)\n")
