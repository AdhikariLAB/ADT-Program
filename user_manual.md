Python setuptools implementation has a lot of options available to properly install a package/library with all necessary
Check 
`python setup.py -h`
`python setup.py --help-commands` to fine tune the installation

### Build:
Instead of installing the package directly one can build the package first to get more (proper adjective ???) control over the installation process. check 
`python setup.py --help-commands` to know more about differnt build options. Also one can build different distributions of the package to quickly install the package in system with similar relavant software configurations and system architecture, using different python package managers. Most usuful ones of them (and also recommended) are the python pip installable pre-built distribution wheel and source distribuition, and can be built by  
`python setup.py bdist_wheel` and  
`python setup.py sdist` respectively and Install them through pip by  
`pip install adt-1.0.0-*.whl`  
`pip install adt-1.0.0-*.tar.gz`  
or you can directly install the packge with the optional requirement i.e h5  
`pip install adt-1.0.0-*.whl[h5]`  
`pip install adt-1.0.0-*.tar.gz[h5]`  

### Install anywhere :  
This package can be installed anywhere in the system. To install it at  an arbitrary directory say DIR, execute  
`python setup.py install --prefix=DIR`  
REMEMBER: DIR must be the absolute path to that dirctory.  
But to work it properly the directory DIR has to be on PYTHONPATH or supports .pth files.  


You can set up the installation directory to support ".pth" files by using one of the approaches described here:
  https://setuptools.readthedocs.io/en/latest/easy_install.html#custom-installation-locations

But the easier option is just to add the directory at the python path.  
set or export the PYTHONPATH using  
`export PYTHONPATH=$PYTHONPATH:DIR/lib/python2.7/site-packages`  
or modify the setup.py by adding this line ( not sure have to test it) anywhere before the setup function  

```python
import sys
sys.path.append('DIR/lib/python2.7/site-packages')` 
```
or the top of the setup.py script add this

```python
import os
os.environ['PYTHONPATH'] += ':'+'DIR/lib/python2.7/site-packages'
```
Now run the installation command to install it in DIR.  
Now as DIR is an arbitrary directory, theres a high chance that DIR is not in your environment path, in that case to use the command line utility `adt` you have to export the location of the excutable i.e.  
`export PATH=$PATH:DIR/bin`  
or set an alias to your `~/.bashrc` file  
`alias = DIR/bin/adt`


### Uninstallation: 
If you have installed the package through pip, then it can be easily uninstalled by  
`pip uninstall adt`.  
On the other hand if you have installed the adt package directly using the setup.py script, then unfortunately, there is no straightforward way to properly remove the package, and you have to manually delete all the files added to your system during the installation process. One way to keep track of the installed files is to provide the  `--record` keyword i.e.  
`python setup.py install --record files.txt`  
now the 'files.txt' will keep the list of all the files added. Now delete all the listed files there to remove the package.


* your system has to have numpy to use this package, though the installer script will try to install numpy from internet, but if you don't have internet connection, the installer will fail. In that case you are requested to install numpy from source(https://docs.scipy.org/doc/numpy/user/building.html) or using other method before hand.

* Installing h5py including the hdf5 api can be preety complicated, thats why h5py installation is optional during the package installation . Install it with any package manager like pip `pip install h5py` or build it from source http://docs.h5py.org/en/latest/build.html 

