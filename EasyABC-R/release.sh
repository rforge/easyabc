#!/bin/sh

DATE=$(date +%Y-%m-%d)
[ $# -eq 1 ] && VERSION=$1 ||
VERSION=$(grep 'Version:'  EasyABC/DESCRIPTION | cut -d' ' -f2)

echo "The package will be build with the version ${VERSION}"
read -p "Are you sure? (y/n) " REPLY
[ $REPLY"" != "y" ]] && {
  read -p "Give the wanted version: " VERSION
}

sed -i -e 's/^Date:.*$/Date: '${DATE}'/' EasyABC/DESCRIPTION EasyABC/man/EasyABC-package.Rd
sed -i -e 's/^Version:.*$/Version: '${VERSION}'/' EasyABC/DESCRIPTION EasyABC/man/EasyABC-package.Rd

R CMD check EasyABC && R CMD build EasyABC
