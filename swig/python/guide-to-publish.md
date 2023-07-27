# How to publish to PyPi

1. Tag and build release on github
Open https://github.com/usgs-coupled/phreeqcrm/releases and press 'Draft a new release'
Press 'Choose a tag' and enter the tag name (ie v0.0.2).  Press the 'Create new tag: v0.0.2 on publish' button.
Enter the 'Release title' (ie PyPi Test Release 0.0.2). Enter Release notes in the 'Describe this release' textbox.
Select 'Set as a pre-release' checkbox. Press 'Publish release' button.

2. When CD workflow finishes download source and wheels artifacts and unzip into dist directory.

3. If not done yet, install build and twine via
```
python -m pip install build twine
```
4. Upload the new files:

TestPyPI (https://test.pypi.org/)
```
python -m twine upload --repository testpypi dist/*
```

PyPI (https://pypi.org/)
```
python -m twine upload dist/*
```