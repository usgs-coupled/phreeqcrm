# Developing PhreeqcRM

## Using a conda environment

### Clone phreeqcrm
```
git clone https://github.com/usgs-coupled/phreeqcrm.git
```

### Create and activate phreeqcrm_dev environment
```
cd phreeqcrm
conda env create -f environment.yml
conda activate phreeqcrm_dev
```

### Configure using CMake Preset
```
cmake --preset phreeqcrm_dev_linux
```

### Build using CMake Preset
```
cmake --build --preset phreeqcrm_dev_linux --parallel $(nproc)
```

### Test using CTest Preset
```
ctest --preset phreeqcrm_dev_linux
```

### Install using CMake Preset
```
cmake --build --preset phreeqcrm_dev_linux --target install
```

## Removing the phreeqcrm_dev environment
```
conda remove --name phreeqcrm_dev --all
```
