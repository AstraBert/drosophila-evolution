command -v conda >/dev/null 2>&1 || {
    echo >&2 "The installation pipeline requires Anaconda/Miniconda but it is not installed. Please check here: https://anaconda.org/ for more details. Aborting."
    exit 1
}

eval "$(conda shell.bash hook)"

for envfile in conda_envs/*.yaml
do
    conda env create -f $envfile
done
