module load roslin/R/4.4.0
module load python/3.11.4
python -m ensurepip --upgrade
python -m venv rnaseq
source rnaseq/bin/activate
pip install -r requirements.txt
