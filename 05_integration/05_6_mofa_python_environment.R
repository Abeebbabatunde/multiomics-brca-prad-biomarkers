## 05_6  MOFA2 Python back-end used for all reported MOFA models
## R MOFA2 1.22.0 + mofapy2 0.7.5 on Python 3.11.9 (reticulate), run_mofa(..., use_basilisk = FALSE).
## basilisk (use_basilisk = TRUE) failed: pyenv_bootstrap_windows() needed git, and git could not resolve github.com.
## Installation (Windows shell):  py -3.11 -m pip install mofapy2 --timeout 300 --retries 10
## RETICULATE_PYTHON must be set BEFORE reticulate is loaded; otherwise reticulate switched to a uv Python 3.12
## without mofapy2.
Sys.setenv(RETICULATE_PYTHON = "C:/Users/USER/AppData/Local/Programs/Python/Python311/python.exe")
library(reticulate)
py_config()
py_module_available("mofapy2")   # must be TRUE
## If a previous run failed: unlink(file.path(od, "MOFA_model.hdf5")) before rerunning
## ("File already has a shared lock" otherwise).
