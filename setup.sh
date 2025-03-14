module load igmm/apps/R/4.4.0
igmm/apps/python/3.12.3
python -m ensurepip --upgrade
python -m venv rnaseq
source rnaseq/bin/activate
pip install -r requirements.txt
