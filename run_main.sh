#!/bin/bash/
cd run/
make clean
make
./vis 3 360 5000000
cd ..

cd 1
make clean
make
(./exampleN01 &)

cd ../2
make clean
make
(./exampleN01 &)

cd ../3
make clean
make
(./exampleN01 &)

cd ../4
make clean
make
(./exampleN01 &)

cd ../5
make clean
make
(./exampleN01 &)

cd ../6
make clean
make
(./exampleN01 &)

cd ../7
make clean
make
(./exampleN01 &)

cd ../8
make clean
make
(./exampleN01 &)

cd ../9
make clean
make
(./exampleN01 &)

cd ../10
make clean
make
(./exampleN01 &)

cd ../11
make clean
make
(./exampleN01 &)

cd ../12
make clean
make
(./exampleN01 &)
cd ..
