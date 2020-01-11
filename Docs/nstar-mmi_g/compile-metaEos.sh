cp makefile/makefile-metaEos sources/makefile
cd sources
make
make clean
cp metaEos.e ../
cd ..
chmod u+x metaEos.e
