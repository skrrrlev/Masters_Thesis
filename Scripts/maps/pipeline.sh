#!/bin/bash

die () {
    echo >&2 "$@"
    exit 9
}

[ "$#" -eq 1 ] || die "1 argument required, $# provided"
[ -d "$1" ] || die "Directory $1 does not exist"

python Scripts/maps/initialize.py $1
python Scripts/maps/dead_pixels.py $1
python Scripts/maps/segmentation_map.py $1
python Scripts/maps/sigma.py $1
python Scripts/maps/crop.py $1
python Scripts/maps/mask.py $1
python Scripts/maps/psf.py $1
python Scripts/maps/galfit_input.py $1