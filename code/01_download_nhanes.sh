#!/usr/bin/env bash
# NHANES bulk downloader: cycles C(2003) D(2005) F(2009) G(2011) H(2013) I(2015) J(2017) L(2021)
RAW="${AGE22_ROOT:-.}/data/raw"
META="${AGE22_ROOT:-.}/data/meta"
mkdir -p "$RAW" "$META"
declare -A YR=( [C]=2003 [D]=2005 [F]=2009 [G]=2011 [H]=2013 [I]=2015 [J]=2017 [L]=2021 )

# candidate file stems per suffix (only suffix matters for component names in modern cycles)
COMMON="DEMO WHQ BMX MCQ RHQ ALQ SMQ PAQ INQ HIQ DBQ"
BIO="BIOPRO CBC CRP HSCRP COT GHB GLU INS TRIGLY HDL TCHOL"
HORM="TST SSTST"
PFAS="PFC PFAS"

download () { # $1=file stem with suffix
  f="$1"
  s="${f##*_}"
  y="${YR[$s]}"
  [ -z "$y" ] && return
  out="$RAW/$f.xpt"
  if [ -s "$out" ]; then echo "SKIP $f"; return; fi
  for attempt in 1 2 3; do
    code=$(curl -s -L --max-time 300 -o "$out.part" -w "%{http_code}" "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/$y/DataFiles/$f.xpt")
    if [ "$code" = "200" ] && [ -s "$out.part" ]; then mv "$out.part" "$out"; echo "OK   $f $(wc -c <"$out")"; return; fi
    sleep 2
  done
  rm -f "$out.part"; echo "MISS $f ($code)"
}

for s in C D F G H I J L; do
  for st in $COMMON; do download "${st}_${s}"; done
  for st in $BIO; do download "${st}_${s}"; done
  for st in $HORM; do
    case "$s" in C|G|H|I|J) download "${st}_${s}";; esac
  done
  for st in $PFAS; do download "${st}_${s}"; done
done

# 2003-2004 individual-level persistent organic pollutants
for f in L28NPB_C L28OCP_C L28DFP_C L24PFC_C SSCHL_C; do download "$f"; done

echo "=== DONE ==="
ls -la "$RAW" | wc -l
