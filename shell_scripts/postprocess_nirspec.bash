# Where the postprocessing scripts are stored. Needs to be set by user.
ROOT=~/Repositories/pdrs4all

set -e
# point the script to where the pipeline output is
RUNDIR="$1"

JWSTVERSION="$(python -c 'import jwst; print(jwst.__version__)')"
PMAPVERSION="$(python -c 'import crds; print(crds.get_context_name("jwst"))')"

NEWDIR=vX.X_auto_"$PMAPVERSION"_"$JWSTVERSION"
echo mkdir -p "$NEWDIR"
mkdir -p "$NEWDIR"

for d in log crf templates
do mkdir -p "$NEWDIR"/"$d"
done

# log files
cp "$RUNDIR"/log/* "$NEWDIR"/log

# crf files (for nirspec, nscleaned output is in separate directory
cp "$RUNDIR"/science/stage3/*crf.fits "$NEWDIR"/crf

# we are now done copying stuff, so cd into the data release directory
cd "$NEWDIR"

# Build the cubes here (separate from main pipeline) because of memory issues. It's also easier
# to customize, and doing it this way we have all our cube builds in one place.

# First, create the association files based on a file pattern (works only for orion data). If
# you use this for other data, then you will have to figure out the appropriate files to build
# each segment.
create_association F100LP_crf crf/jw01288003001_0312?_*crf.fits
create_association F170LP_crf crf/jw01288003001_0311?_*crf.fits
create_association F290LP_crf crf/jw01288003001_0310?_*crf.fits

# Default resolution cubes (but with custom PA to optimize cube size)
mkdir -p cubes/default
for ASN in F???LP_crf_asn.json
do strun cube_build $ASN --output_dir cubes/default --cube_pa=250.42338204969806
done

# Cubes matched to MIRI MRS Ch1
mkdir -p cubes/ch1wcs
for ASN in F???LP_crf_asn.json
do strun cube_build $ASN --output_dir cubes/ch1wcs --cube_pa=250.42338204969806 --scalexy 0.13 --ra_center 83.83535169747073 --dec_center -5.419729392828107  --nspax_x 213 --nspax_y 41
done

# WCS correction needs to happen here
python3 -m pdrs4all.postprocess.nirspec_wcs_calibrate_stitch cubes/default/*s3d.fits --output_dir cubes/default_wcscorr

# templates from both cubes, both stitched and not stitched
extract_templates "$ROOT"/regions/aper_T_DF_extraction.reg cubes/default/*s3d.fits --template_names "HII" "Atomic" "DF3" "DF2" "DF1" -o templates/default_nostitch.ecsv

extract_templates "$ROOT"/regions/aper_T_DF_extraction.reg cubes/default/*s3d.fits --template_names "HII" "Atomic" "DF3" "DF2" "DF1" --apply_offsets --reference_segment 2 -o templates/default_addstitch.ecsv

extract_templates "$ROOT"/regions/aper_T_DF_extraction.reg cubes/default_wcscorr/Level3* --template_names "HII" "Atomic" "DF3" "DF2" "DF1" -o templates/default_wcscorr_nostitch.ecsv

extract_templates "$ROOT"/regions/aper_T_DF_extraction.reg cubes/default_wcscorr/Level3* --template_names "HII" "Atomic" "DF3" "DF2" "DF1" --apply_offsets --reference_segment 2 -o templates/default_wcscorr_addstitch.ecsv
