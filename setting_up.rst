Running ColRadPy
================

Go to desired directory and in the terminal, run

``git clone https://github.com/johnson-c/ColRadPy.git``

and in ColRadPy directory,

``git init``

to initialize git repo. Use ``git branch -a`` in ColRadPy to make sure you have the HEAD
file.

Then, also in the ColRadPy directory,

``python3 setup.py build``

``python3 setup.py install``

Finally, when you append your path in your python file like in the tutorial, use
``sys.path.append('path')`` , where ``'path'`` is the path to where the colradpy folder is (in the
ColRadPy folder if you haven’t changed anything).
