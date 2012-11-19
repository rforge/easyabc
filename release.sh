#!/bin/sh

DATE=$(date +%Y-%m-%d)
[ $# -eq 1 ] && VERSION=$1 ||
VERSION=$(grep 'Version:'  pkg/DESCRIPTION | cut -d' ' -f2)

echo "The package will be build with the version ${VERSION}"
read -p "Are you sure? (y/n) " REPLY
[ $REPLY"" != "y" ]] && {
  read -p "Give the wanted version: " VERSION
}

sed -i -e 's/^Date:.*$/Date: '${DATE}'/' pkg/DESCRIPTION pkg/man/EasyABC-package.Rd
sed -i -e 's/^Version:.*$/Version: '${VERSION}'/' pkg/DESCRIPTION pkg/man/EasyABC-package.Rd
sed -i -e 's/\(\\date{\\texttt{EasyABC} version \)[^,]*,/\1 '${VERSION}'/' vignettes/EasyABC.Rnw

./updateVignette.sh && \
  R CMD build pkg && \
  R CMD check --as-cran EasyABC_${VERSION}.tar.gz
