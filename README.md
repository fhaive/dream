# DREAM  - Druggability Evaluation of Human Complex Diseases
This is the source code of "DREAM: an R package for druggability evaluation of human complex diseases".
DREAM focuses on assessing the druggability of specific pathological conditions by leveraging robust drug repurposing predictions.
It provides a powerful tool for identifying potential sets of drugs putatively suitable for combination therapy.
## Environment
To create the necessary environment use the following command 
> :warning: note: the environment has only been tested on Linux! For other systems please refer to the Docker section
```
conda env create --file environment.yml -n dream
```
### Other dependencies
After creating the environment run the following commands to complete the installation of the required dependencies

```
conda run -n dream R -e "install.packages('TopKLists', dependencies=FALSE, repos='http://cran.us.r-project.org')"
conda run -n dream R -e "devtools::install_github('filosi/nettools', dependencies=FALSE, upgrade_dependencies=FALSE)"
conda run -n dream R -e "install.packages('IRkernel', dependencies=FALSE, repos='http://cran.us.r-project.org'); IRkernel::installspec()"
conda run -n dream R -e "install.packages('http://localhost:9000/DREAM_0.1.0.tar.gz', repos=NULL)"
conda run -n dream pip install http://localhost:9000/dist/dream-0.1.dev0-py3-none-any.whl
```
## Use the environment
The environment provides a Jupyter server in which it is possible to run R code and perform analyses with DREAM.
To start the Jupyter server run the following command:
```
conda run -n dream jupyter lab --NotebookApp.token='' --NotebookApp.password=''
```
And browse to the address `http://localhost:8888/lab`
It is possible to either start from a new notebook, or open the provided notebooks that shows how to use the functionalities provided by the package.
Specifically, the notebook `dream_vignette.ipynb` shows how the analysis reported as the case study described in Supplementary File 2 of the publication was performed.

## Docker
We also provide a ready made docker environment to try the package that can be launched using the command
```
docker run -it --rm -p8888:8888 -e UID=$(id -u) -e GID=$(id -g) -v$(pwd):/workspace fhaive/dream:latest
```
This will start the same jupyter environment within the docker that can be accessed in the same manner browsing the page `http://localhost:8888/lab`
In addition, the current directory in which the command is executed will be mapped inside the docker container, facilitating data exchange.

### Build the docker image
To build the docker image from scratch, use the following commands:
```
# clone the repository
git clone git@github.com:fhaive/dream.git
cd dream/
# follow the instructions above to create the environment
# activate the environment
conda activate dream
# build the python package
cd Python_package/
python setup.py bdist_wheel
# build the R package
cd ../R_package/
R CMD build .
# build the docker
cd ..
docker build -t dream .
# after building, run the container
docker run -it --rm -p8888:8888 -e UID=$(id -u) -e GID=$(id -g) -v$(pwd):/workspace dream
```