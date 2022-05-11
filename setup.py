from setuptools import setup


setup(
    name = "MTLib",
    author = "Ditlev Frickmann @skrrrlev",
    author_email = "reimer.frickmann@gmail.com",
    description = "Master Thesis Library",
    version = "0.0.0",
    license = "MIT",
    url = "https://github.com/skrrrlev/Master_Thesis/MTLib",  
    packages=[
        'MTLib',
        'MTLib.data',
        'MTLib.data.stardust',
        'MTLib.data.cosmos',
        'MTLib.data.dataclasses',
        'MTLib.cosmos',
        'MTLib.plots',
        'MTLib.sfr',
        'MTLib.fitstools',
        'MTLib.files',
    ],
    package_data={
        'MTLib.cosmos':['*.ini']
    },
    include_package_data=True,
    install_requires = [
        'numpy==1.19.5',
        'astropy==4.1',
        'C4S@git+https://github.com/skrrrlev/Cataloguer-4-Stardust.git'
    ],
)