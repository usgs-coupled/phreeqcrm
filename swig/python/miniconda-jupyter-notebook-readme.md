## Setup miniconda to use phreeqcrm in a jupyter notebook (~5 GB)

Install miniconda (https://docs.conda.io/en/latest/miniconda.html)

![Welcome](images/miniconda-0.png)

![License](images/miniconda-1.png)

![Just Me](images/miniconda-2.png)

![Destination](images/miniconda-3.png)

![Options](images/miniconda-4.png)

![Finish](images/miniconda-5.png)

### Open 'Anaconda Powershell Prompt (miniconda3)' from the 'Miniconda3 (64-bit)' menu item

#### Create phreeqcrm environment:
```
conda create --name phreeqcrm python=3.11 --yes
```

#### Activate phreeqcrm environment
```
conda activate phreeqcrm
```

#### Install prerequisites (this will become unnecessary when on https://pypi.org)
```
conda install numpy PyYAML --yes
```

#### Install matplotlib (optional but recommended)
```
conda install matplotlib --yes
```

#### Install phreeqcrm
```
pip install -i https://test.pypi.org/simple/ phreeqcrm
```

#### Install Jupyter Notebook
```
conda install -c conda-forge notebook --yes
```

#### Change to ex11-advect directory
```
cd .\swig\python\ex11-advect\
```

#### Start jupyter notebooks
```
jupyter notebook
```

#### In the web browser double-click the `ex11-advect.ipynb` file

#### To test the install use the menu item:
```
Kernel -> Restart Kernel and Run All Cells... -> Restart
```

## Uninstall miniconda
(see https://stackoverflow.com/questions/29596350/how-to-uninstall-mini-conda-python)

### Open 'Anaconda Powershell Prompt (miniconda3)' from the 'Miniconda3 (64-bit)' menu item

#### Activate phreeqcrm environment
```
conda activate phreeqcrm
```

#### Install anaconda-clean
```
conda install anaconda-clean
```

#### Run anaconda-clean
```
anaconda-clean --yes
```

#### Uninstall Miniconda

Open Control Panel -> Programs -> Uninstall a program and select 'miniconda``
