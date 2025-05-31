# monitoring_stuff

## 📌 Overview

This repository contains a modular workflow for downloading, processing, and analyzing seismic data using ambient noise techniques. The project is designed to support the full monitoring pipeline: from raw data acquisition to estimation of medium changes.

Scripts are organized into logically ordered folders for clarity and scalability. Parameters such as paths, filter settings, and component choices are centrally configured via `config.py`.

---

## 🗂 Project Structure

```
monitoring_stuff/
├── A_download/          # Download waveform data from FDSN services
├── B_preprocess/        # Clean, filter, and quality control of traces
├── C_correlations/      # Compute auto- and cross-correlations between station pairs
├── D_stacks/            # Temporal stacking of correlations (e.g. daily, weekly)
├── E_changes/           # Quantify changes over time
├── files/               # Input metadata (station lists, coordinates, etc.)
├── config.py            # Central configuration for all modules
└── README.md            # Project documentation
```

---

## ⚙ Configuration

All parameters (input/output paths, filters, sampling rates, etc.) are managed in `config.py`. Adjust this file before running any module to match your local directory structure and preferences.

---

## 🛠 Requirements

- Python ≥ 3.7
- ObsPy
- NumPy
- Pandas

---

## 🧪 Output Notes

Due to the potentially large size of waveform and correlation data, outputs (e.g. SAC files) are **not stored inside the repository**. Output paths are user-defined in `config.py` and can point to external storage locations.

---

## 🚧 Development Status

Currently, the download module is functional and ready to use. Additional modules for preprocessing, correlation, and medium changes analysis will be added progressively.

---

## ✉ Contact

For questions, suggestions, or contributions, please contact:
Pilar Sánchez^2-Pastor
psanchez@geo3bcn.csic.es
