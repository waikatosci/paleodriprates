#!/usr/bin/env bash
# fetch_external.sh — retrieve the two third-party gridded inputs Fig 2 needs.
#
# These are published products, not our data: they are cited, not archived. Fetching
# beats vendoring 152 MB of someone else's raster into a repo.
# Compare the checksums afterwards against external/README.md.
set -euo pipefail
cd "$(dirname "$0")"

echo "1/2  GPCP v2.3 monthly precipitation (NOAA PSL, ~20 MB) ..."
curl -fL -o precip.mon.mean.nc \
  https://psl.noaa.gov/thredds/fileServer/Datasets/gpcp/precip.mon.mean.nc

echo "2/2  Natural Earth 10m shaded relief (~131 MB unzipped) ..."
curl -fL -o SR_LR.zip https://naturalearth.s3.amazonaws.com/10m_raster/SR_LR.zip
unzip -o SR_LR.zip 'SR_LR/SR_LR.tif' && mv SR_LR/SR_LR.tif . && rm -rf SR_LR SR_LR.zip

echo
echo "Verify against the sha256 values in README.md:"
sha256sum precip.mon.mean.nc SR_LR.tif

echo
echo "NOTE: Fig 2 also fetches Natural Earth COASTLINES at render time via cartopy,"
echo "      so it needs network access when you run it — not just now."
