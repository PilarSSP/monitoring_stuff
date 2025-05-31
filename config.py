# config.py

CONFIG = {
    # Processing and filter settings for removing the instrument response
    'bandpass_filter': (0.01, 0.05, 10, 20),
    'resample_rate': 20,

    # Directory structure
    'output_base_dir': '/Users/pilarsanchez/Desktop/tests',
    'window_hours': 1,
    'window_label': '1h',

    # Input files (located inside: files/)
    'station_list_file': 'in_stations.txt',
    'station_info_file': 'stations_info.txt',

    # Channel selection and priority order if multiple sensors are available for the same station
    'component': 'Z',
    'channel_priority': ['HHZ', 'BHZ', 'LHZ'],

    # Data source
    'fdsn_url': 'http://service.iris.edu/',
    
    # Maximum number of parallel download threads. Adjust according to your CPU cores
    'max_workers': 4
}
