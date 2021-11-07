#!/bin/sh
xmake g --network=private
xrepo import -i ./packages
xmake -y
