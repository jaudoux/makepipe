#! /bin/bash

VERSION=$(grep -oP "version\s+\K(\S+)" makepipe)
tar -cvzf makepipe-$VERSION.tar.gz bricks makepipe
