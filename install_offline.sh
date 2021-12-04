#!/bin/sh

# Check installation of Xmake
if [ $(which xmake) ];
then
rm -rf ./xmake
rm -rf ./build-artifacts
rm -rf ./xmake-repo
else
cd ./xmake
make build 
./scripts/get.sh __local__ __install_only__
source ~/.xmake/profile
cd ..
rm -rf ./xmake
mv ./build-artifacts ~/.xmake/repositories/
mv ./xmake-repo ~/.xmake/repositories/
fi

if [ $(which MMC) ];
then 
# Uninstall old version
xmake uninstall --installdir=$HOME/.local/MMC
rm $HOME/.local/bin/MMC
else 
# Install packages
xmake g --network=private
xrepo import -i ./packages cmake scnlib fmt ctre spdlog toml++
fi

# Install MMC
xmake -y
xmake install -o $HOME/.local/MMC
ln -s $HOME/.local/MMC/bin/MMC $HOME/.local/bin/MMC