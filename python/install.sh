nppath=`python -c 'import numpy; print(numpy.get_include())'`
export PYTHONPATH=$PYTHONPATH:./lib/python/
CFLAGS="-I $nppath/numpy" python setup.py install --home=./
python setup.py install
cp build/lib.*/dimension/hidalgo/_gibbs* ./dimension/_gibbs.so

