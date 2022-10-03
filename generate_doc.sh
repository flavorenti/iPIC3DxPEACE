#!/bin/bash

cd ..

doxygen Doxyfile_iPIC3DxPEACE

cd ./flavorenti.github.io

git add .

git commit -m 'ciao'

git push
