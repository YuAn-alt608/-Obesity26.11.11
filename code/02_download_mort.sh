#!/usr/bin/env bash
RAW="${AGE22_ROOT:-.}/data/raw"
mkdir -p "$RAW/mort"
for c in 2003_2004 2005_2006 2009_2010 2011_2012 2013_2014 2015_2016 2017_2018; do
  f="NHANES_${c}_MORT_2019_PUBLIC.dat"
  out="$RAW/mort/$f"
  [ -s "$out" ] && { echo "SKIP $f"; continue; }
  code=$(curl -s -L --max-time 300 -o "$out" -w "%{http_code}" "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/$f")
  echo "$code $f $(wc -c <"$out")"
done
ls -la "$RAW/mort"
