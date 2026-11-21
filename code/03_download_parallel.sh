#!/usr/bin/env bash
RAW="${AGE22_ROOT:-.}/data/raw"
mkdir -p "$RAW" "$RAW/geo" "${AGE22_ROOT:-.}/data/meta"
LIST="${AGE22_ROOT:-.}/data/meta/download_manifest.txt"
: > "$LIST"

add () { echo "$1|$2" >> "$LIST"; }   # file|year

# ---------- 2003-2004 (cycle C): legacy component filenames ----------
for f in DEMO_C WHQ_C BMX_C MCQ_C RHQ_C ALQ_C SMQ_C PAQ_C HIQ_C DBQ_C \
         L40_C L25_C L11_C L10_C L10AM_C L13_C L13AM_C L06COT_C SMQMEC_C \
         L28NPB_C L28OCP_C L28DFP_C L24PFC_C SSCHL_C RXQ_RX_C; do add "$f" 2003; done

# ---------- 2005-2006 (D) ----------
for f in DEMO_D WHQ_D BMX_D MCQ_D RHQ_D ALQ_D SMQ_D PAQ_D HIQ_D DBQ_D INQ_D \
         BIOPRO_D CBC_D CRP_D GHB_D GLU_D TCHOL_D HDL_D TRIGLY_D COT_D PFC_D RXQ_RX_D; do add "$f" 2005; done

# ---------- 2009-2010 (F) ----------
for f in DEMO_F WHQ_F BMX_F MCQ_F RHQ_F ALQ_F SMQ_F PAQ_F HIQ_F DBQ_F INQ_F \
         BIOPRO_F CBC_F CRP_F GHB_F GLU_F INS_F TCHOL_F HDL_F TRIGLY_F COT_F PFC_F RXQ_RX_F; do add "$f" 2009; done

# ---------- 2011-2012 (G) ----------
for f in DEMO_G WHQ_G BMX_G MCQ_G RHQ_G ALQ_G SMQ_G PAQ_G HIQ_G DBQ_G INQ_G \
         BIOPRO_G CBC_G CRP_G GHB_G GLU_G INS_G TCHOL_G HDL_G TRIGLY_G COTNAL_G PFC_G TST_G RXQ_RX_G; do add "$f" 2011; done

# ---------- 2013-2014 (H) ----------
for f in DEMO_H WHQ_H BMX_H MCQ_H RHQ_H ALQ_H SMQ_H PAQ_H HIQ_H DBQ_H INQ_H \
         BIOPRO_H CBC_H CRP_H HSCRP_H GHB_H GLU_H INS_H TCHOL_H HDL_H TRIGLY_H COT_H PFAS_H TST_H RXQ_RX_H; do add "$f" 2013; done

# ---------- 2015-2016 (I) ----------
for f in DEMO_I WHQ_I BMX_I MCQ_I RHQ_I ALQ_I SMQ_I PAQ_I HIQ_I DBQ_I INQ_I \
         BIOPRO_I CBC_I HSCRP_I GHB_I GLU_I INS_I TCHOL_I HDL_I TRIGLY_I COT_I PFAS_I TST_I RXQ_RX_I; do add "$f" 2015; done

# ---------- 2017-2018 (J) ----------
for f in DEMO_J WHQ_J BMX_J MCQ_J RHQ_J ALQ_J SMQ_J PAQ_J HIQ_J DBQ_J INQ_J \
         BIOPRO_J CBC_J HSCRP_J GHB_J GLU_J INS_J TCHOL_J HDL_J TRIGLY_J COT_J PFAS_J SSTST_J RXQ_RX_J; do add "$f" 2017; done

# ---------- 2021-2023 (L) ----------
for f in DEMO_L WHQ_L BMX_L MCQ_L RHQ_L ALQ_L SMQ_L PAQ_L HIQ_L DBQ_L INQ_L \
         BIOPRO_L CBC_L HSCRP_L GHB_L GLU_L INS_L TCHOL_L HDL_L TRIGLY_L COT_L PFAS_L RXQ_RX_L; do add "$f" 2021; done

fetch () {
  line="$1"; f="${line%%|*}"; y="${line##*|}"
  out="$RAW/$f.xpt"
  [ -s "$out" ] && { echo "SKIP $f"; return; }
  url="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/$y/DataFiles/$f.xpt"
  code=$(curl -s -f -L --max-time 240 --retry 2 -o "$out.part" -w "%{http_code}" "$url")
  if [ "$code" = "200" ] && [ -s "$out.part" ]; then mv "$out.part" "$out"; echo "OK   $f $(wc -c <"$out")";
  else rm -f "$out.part"; echo "MISS $f $code"; fi
}
export -f fetch
export RAW

cat "$LIST" | xargs -P 8 -I{} bash -c 'fetch "$@"' _ {} > "${AGE22_ROOT:-.}/logs/03_dl.log" 2>&1
echo "=== downloaded: $(ls $RAW/*.xpt 2>/dev/null | wc -l) files ==="
