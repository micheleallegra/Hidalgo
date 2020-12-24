nppath=`python -c 'import numpy; print(numpy.get_include())'`
export PYTHONPATH=$PYTHONPATH:./lib/python/
CFLAGS="-I $nppath/numpy" python setup.py install --home=./
python setup.py build
cp build/lib.macosx-10.9-x86_64-3.8/_gibbs.cpython-38-darwin.so ./dimension/_gibbs.so
