# Setup

In this page, I want to explain the process that I've taken to setup the environment for this project.

## Cloning
If you want all of the scripts and data please clone this entire repository to your machine. Do note that I wrote this on a Mac, however, it should also work for Windows and Linux. 

If you're unfamiliar with git, or simply really annoyed with the console-based approach, I would recomend installing a git giu. I really like [SmartGit](https://www.syntevo.com/smartgit/), however since it is not free, I usually use [GitKraken](https://www.gitkraken.com/).

## Python installation
This project uses the [Stardust SED fitting tool](https://github.com/VasilyKokorev/stardust), which is a python based approach to fitting SEDs of galaxies.

Unfortunately _Stardust_ was developed for Python version 3.6.10. If you have a Python version >3.6.ZZ, _Stardust_ does not work.
You can check your version in the terminal:
```console
python -V
>> Python X.Y.ZZ
```

If you're version is >3.6.ZZ then you might want to consider installing _pyenv_. It is a package that can easily enable you to have multiple installations of python on your computer. I recommend to go through [this tutorial](https://k0nze.dev/posts/install-pyenv-venv-vscode/); It has sections for every operating system.

When you have _pyenv_ installed, open a terminal and install Python version 3.6.15. (Only the latest version of older packages are supported, but _Stardust_ works fine with this version.)
```console
pyenv install 3.6.15
```

After it has installed, go to the folder where you cloned this git repository and enable Python 3.6.15 for that folder:
```console
pyenv local 3.6.15
```

When you invoke python from the command line, you should now see that it defaults to the 3.6.15 version. 
```console
python -V
>> Python 3.6.15
```

## Virtual environment
Often when working with something in Python, it is beneficial to setup a virtual environment. It is basically a clean sheet of python, where only the built-in modules are installed. Other packages that might come in handy, can easily be installed with _pip_.

This can simply be done by invoking the _virtualenv_ python package from the console. 
The package is usually installed already
```console
python -m pip install virtualenv
```

It is customary to name your virtual environment _venv_, but you can choose whatever name you prefer. Here, I create a virtual environment called _venv_ in the current directory (the root directory of the project), using the console.
Then, initialise a virtual environment with the name _venv_:
```console
python -m virtualenv venv
```

Now, activate the environment
```console
source venv/bin/activate    # (mac/linux)
call venv\Scripts\activate  # (windows)
```

If you're in an IDE, you should also select the python executable path to
```console
venv/bin/python    # (mac/linux)
venv\Scripts\python  # (windows)
```

Now, when you invoke python, it will invoke the executable in the virtual environment. 

## Requirements
I set up a list of requirements in the root folder of this project. These are packages used by the library and the scripts. Install them to your virtual environment using the following console command from the root of the project, while having the virtual enviroment activated,
```console
python -m pip install -r requirements.txt
```
This list also includes the _Stardust_ package, which will be installed directly from its GitHub page. Moreover, it also includes the package [Cataloguer 4 Stardust (C4S)](https://github.com/skrrrlev/Cataloguer-4-Stardust) that I created to make it easier to work with _Stardust_.

## MTLib: Master Thesis Library
I've created a library of methods and classes used in the scripts of this project. Install them to your virtual environment using the _setup.py_ script in the root folder:
```console
python -m pip install .
```