cp makefile/makefile-crustSna sources/makefile
cd sources
make
make clean
cp crustSna.e ../
cd ..
chmod u+x crustSna.e
