#!/usr/bin/env bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

rm -rf /tmp/easistrain_results

python -m easistrain.EDD.calibrationEDD $SCRIPT_DIR/calibEDD.yml

python -m easistrain.EDD.angleCalibEDD $SCRIPT_DIR/angleCalibEDD.yml

python -m easistrain.EDD.fitEDD $SCRIPT_DIR/fitEDD_OR1.yml

python -m easistrain.EDD.fitEDD $SCRIPT_DIR/fitEDD_OR2.yml

python -m easistrain.EDD.coordTransformation $SCRIPT_DIR/coordTransform_OR1.yml

python -m easistrain.EDD.coordTransformation $SCRIPT_DIR/coordTransform_OR2.yml

python -m easistrain.EDD.regroupPoints $SCRIPT_DIR/regroupPoints.yml

python -m easistrain.EDD.preStraind0cstEDD $SCRIPT_DIR/preStrain.yml

python -m easistrain.EDD.strainStressd0cstEDD $SCRIPT_DIR/strain.yml
