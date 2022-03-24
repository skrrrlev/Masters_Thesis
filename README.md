# Master Thesis Git Repository
This repository contains all the code I write and data used/produced for my Master's thesis.


## Setup

First, clone the repository to your place of preference, and open a console inside the root of the repository.

Unfortunately _Stardust_ was developed for Python version 3.6.10. If you have a Python version >3.6.ZZ, _Stardust_ does not work. You can check your version in the terminal:
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

Running Python from inside this folder will then invoke this version. Next, to run _Stardust_ you should install it and its dependencies. I have compiled a list in the _requirements.txt_ document in the root folder that can be installed using _pip_. However, I would suggest you set up a virtual environment first:

First, make sure you have the _virtualenv_ python package installed:
```console
python -m pip install virtualenv
```

Then, initialise a virtual environment with the name _venv_:
```console
python -m virtualenv venv
```

Now, activate the environment
```console
source venv/bin/activate    # (mac/linux)
call venv\Scripts\activate  # (windows)
```

Finally, install the requirements:
```console
python -m pip install -r requirements.txt   # installs all the dependecies + Stardust.
python -m pip install .                     # installs the library I built for this project.
```

You are now ready to proceed.