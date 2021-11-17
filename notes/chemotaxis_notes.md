
### Setup:
1) clone GitHub repo to directory of your choice.

```
git clone https://github.com/AndersenLab/chemotaxis-cli.git
```

2) Build virtual environment with python 2.7.11. if not already built.
	
```
pyenv virtualenv 2.7.11 "chemotaxis"
```

3) Nativagte to the repo directory then set the environment.

```

cd <your_directory>
pyenv local "chemotaxis"
eval "$(pyenv init -)" - tim edit 20211116
```

4) install requirements for package

```
pip install -r requirements.txt
```

5) install setup utility

```
python setup.py install
```


### After running script for the first time:
1) put images in the chemotaxis directory then change directory to image directory

2) Run ct script. The code below will run ct script on all .jpg files in the image directory. 
# --radius tells the script to use a radius of  865 pixels (you may change this), -d means debug mode,

Measure diameter of plate or inner circle in imageJ, then divide by two to find radius for cicles.

# --header means make txt file with a header and call it results.txt

`ct *.jpg  --radius 865  -d --header > results.txt`

# notes:
        1) Find radius in imageJ. May need to add additional crop with --crop XX pixels
        2) 




