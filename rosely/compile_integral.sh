swig -python -c++ src/integral.i
cp src/integral.py .
python3 setup.py build_ext
cp build/lib.macosx-10.6-intel-3.6/rosely/_integral.cpython-36m-darwin.so .
