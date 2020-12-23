nppath=`python -c 'import numpy; print(numpy.__path__[0])'`
nppath=$nppath/core/include/numpy
mkdir lib
mkdir lib/python
export PYTHONPATH=$PYTHONPATH:./lib/python/
CFLAGS="-I $nppath" python setup.py install --home=./ 
python setup.py build
gs=`find build/ -name _gibbs.so`
cp $gs dimension/
