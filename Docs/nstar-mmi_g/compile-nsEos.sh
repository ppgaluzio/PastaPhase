cp makefile/makefile-nsEos sources/makefile
cd sources
make
make clean
cp nsEos.e ../
cd ..
chmod u+x nsEos.e
