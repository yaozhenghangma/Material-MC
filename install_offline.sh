#!/bin/sh

if [ $(which MMC) ];
then 
# Uninstall old version
rm $HOME/.local/bin/MMC
fi

# Install MMC
cmake . -DCMAKE_CXX_COMPILER=`which g++` -DCMAKE_C_COMPILER=`which gcc`
cmake --build .
ln build/MMC ~/.local/bin/MMC