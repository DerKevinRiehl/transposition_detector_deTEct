#!/bin/bash
mkdir -p $PREFIX/bin
mkdir -p $PREFIX/share
mkdir -p $PREFIX/share/TranspositionEventDetector_deTEct
cp $RECIPE_DIR/*.py $PREFIX/share/TranspositionEventDetector_deTEct
cp $RECIPE_DIR/transposition_deTEct $PREFIX/bin
chmod +x $PREFIX/bin/transposition_deTEct
