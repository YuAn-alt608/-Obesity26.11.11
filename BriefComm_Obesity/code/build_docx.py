# -*- coding: utf-8 -*-
"""build_docx.py -- assemble the submission files as .docx"""
import os
import re
import csv
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.enum.section import WD_SECTION
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MS = os.path.join(BASE, "manuscript")
SUP = os.path.join(BASE, "supplement")


# ---------------------------------------------------------------- helpers
def add_page_number_footer(section):
    p = section.footer.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run()
    fld = OxmlElement("w:fldSimple")
    fld.set(qn("w:instr"), "PAGE")
    run._r.addnext(fld)


def add_line_numbers(section):
    sectPr = section._sectPr
    ln = OxmlElement("w:lnNumType")
    ln.set(qn("w:countBy"), "1")
    ln.set(qn("w:restart"), "continuous")
    sectPr.append(ln)


def base_style(doc, size=12):
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(size)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    pf = st.paragraph_format
    pf.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    pf.space_before = Pt(0)
    pf.space_after = Pt(0)


SUPERSCRIPT_MAP = {
    "\u2070": "0", "\u00b9": "1", "\u00b2": "2", "\u00b3": "3", "\u2074": "4",
    "\u2075": "5", "\u2076": "6", "\u2077": "7", "\u2078": "8", "\u2079": "9",
    "\u207b": "\u2212", "\u207a": "+",
}


def _emit(p, text, size, bold=False, italic=False, mono=False):
    """Write text, turning Unicode superscripts into real superscript runs."""
    font = "Courier New" if mono else "Times New Roman"
    buf = ""

    def flush():
        nonlocal buf
        if buf:
            r = p.add_run(buf)
            r.bold, r.italic = bold, italic
            r.font.size = Pt(size)
            r.font.name = font
            buf = ""

    for ch in text:
        if ch in SUPERSCRIPT_MAP:
            flush()
            r = p.add_run(SUPERSCRIPT_MAP[ch])
            r.font.superscript = True
            r.bold, r.italic = bold, italic
            r.font.size = Pt(size)
            r.font.name = "Times New Roman"
        else:
            buf += ch
    flush()


def add_runs(p, text, size=12, bold_all=False, italic=False):
    """Handle **bold**, *italic* and `code` inline markers."""
    for tok in re.split(r"(\*\*.+?\*\*|`[^`]+`|\*[^*]+\*)", text):
        if not tok:
            continue
        if tok.startswith("**") and tok.endswith("**") and len(tok) > 4:
            inner = tok[2:-2]
            for sub in re.split(r"(\*[^*]+\*)", inner):
                if not sub:
                    continue
                if sub.startswith("*") and sub.endswith("*") and len(sub) > 2:
                    _emit(p, sub[1:-1], size, bold=True, italic=True)
                else:
                    _emit(p, sub, size, bold=True, italic=italic)
        elif tok.startswith("`") and tok.endswith("`") and len(tok) > 2:
            _emit(p, tok[1:-1], size, bold=bold_all, mono=True)
        elif tok.startswith("*") and tok.endswith("*") and len(tok) > 2:
            _emit(p, tok[1:-1], size, bold=bold_all, italic=True)
        else:
            _emit(p, tok, size, bold=bold_all, italic=italic)
    return p


def page_setup(doc):
    for s in doc.sections:
        s.page_width = Inches(8.5)
        s.page_height = Inches(11)
        for attr in ("top_margin", "bottom_margin", "left_margin", "right_margin"):
            setattr(s, attr, Inches(1))


def md_table_to_docx(doc, rows, font=8.5):
    ncol = len(rows[0])
    t = doc.add_table(rows=len(rows), cols=ncol)
    t.style = "Table Grid"
    t.autofit = True
    for i, row in enumerate(rows):
        for j, cell in enumerate(row):
            c = t.cell(i, j)
            c.text = ""
            p = c.paragraphs[0]
            p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            p.paragraph_format.space_before = Pt(1)
            p.paragraph_format.space_after = Pt(1)
            add_runs(p, cell, size=font, bold_all=(i == 0))
    return t


def parse_table(lines):
    rows = []
    for ln in lines:
        if set(ln.replace("|", "").strip()) <= set("-: "):
            continue
        cells = [c.strip() for c in ln.strip().strip("|").split("|")]
        rows.append(cells)
    return rows


def render_markdown(doc, md_path, body_size=12, title_size=14):
    lines = open(md_path, encoding="utf-8").read().split("\n")
    i = 0
    in_refs = False
    while i < len(lines):
        ln = lines[i].rstrip()
        if not ln.strip():
            i += 1
            continue
        if ln.strip() == "---":
            i += 1
            continue
        if ln.startswith("|"):
            block = []
            while i < len(lines) and lines[i].startswith("|"):
                block.append(lines[i].rstrip())
                i += 1
            md_table_to_docx(doc, parse_table(block))
            continue
        if ln.startswith("# "):
            p = doc.add_paragraph()
            p.paragraph_format.space_after = Pt(12)
            r = p.add_run(ln[2:].strip())
            r.bold = True
            r.font.size = Pt(title_size)
            i += 1
            continue
        if ln.startswith("## "):
            txt = ln[3:].strip()
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(12)
            p.paragraph_format.space_after = Pt(4)
            r = p.add_run(txt)
            r.bold = True
            r.font.size = Pt(body_size)
            in_refs = txt.lower().startswith("reference")
            i += 1
            continue
        if ln.startswith("- "):
            p = doc.add_paragraph(style="List Bullet")
            p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
            add_runs(p, ln[2:], size=body_size)
            i += 1
            continue
        if re.match(r"^\d+\.\s", ln) and in_refs:
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Inches(0.3)
            p.paragraph_format.first_line_indent = Inches(-0.3)
            add_runs(p, ln, size=body_size)
            i += 1
            continue
        # ordinary paragraph
        p = doc.add_paragraph()
        p.paragraph_format.space_after = Pt(0)
        add_runs(p, ln, size=body_size)
        i += 1


# ---------------------------------------------------------------- 1. main
doc = Document()
base_style(doc, 12)
page_setup(doc)
add_page_number_footer(doc.sections[0])
add_line_numbers(doc.sections[0])
render_markdown(doc, os.path.join(MS, "main_manuscript.md"))
out1 = os.path.join(MS, "Main_Manuscript.docx")
doc.core_properties.author = "Zizheng Yuan"
doc.core_properties.title = ("Three artefacts manufacture a \"weight-loss toxin rebound\" "
                             "alarm, and total testosterone masks a fall in the free fraction")
doc.save(out1)

# ---------------------------------------------------------------- 2. supplement
doc = Document()
base_style(doc, 11)
page_setup(doc)
add_page_number_footer(doc.sections[0])
render_markdown(doc, os.path.join(SUP, "supplementary_material.md"),
                body_size=11, title_size=13)
out2 = os.path.join(SUP, "Supplementary_Material.docx")
doc.save(out2)

print("saved:", out1)
print("saved:", out2)

# ---------------------------------------------------------------- 3. cover letter
doc = Document()
st = doc.styles["Normal"]
st.font.name = "Times New Roman"
st.font.size = Pt(12)
st.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
st.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
st.paragraph_format.space_after = Pt(8)
page_setup(doc)

for ln in open(os.path.join(MS, "cover_letter.md"), encoding="utf-8").read().split("\n"):
    ln = ln.rstrip()
    if not ln.strip() or ln.strip() == "---":
        doc.add_paragraph()
        continue
    if ln.startswith("# "):
        p = doc.add_paragraph()
        r = p.add_run(ln[2:].strip()); r.bold = True
        continue
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    add_runs(p, ln, size=12)
out3 = os.path.join(MS, "Cover_Letter.docx")
doc.save(out3)

# ---------------------------------------------------------------- 4. STROBE
doc = Document()
base_style(doc, 10)
page_setup(doc)
add_page_number_footer(doc.sections[0])
p = doc.add_paragraph()
r = p.add_run("STROBE checklist")
r.bold = True; r.font.size = Pt(13)
p = doc.add_paragraph()
add_runs(p, "Manuscript: Three artefacts manufacture a \"weight-loss toxin rebound\" alarm, and total testosterone masks a fall in the free fraction.", size=10)
rows = list(csv.reader(open(os.path.join(SUP, "STROBE_checklist.csv"), encoding="utf-8")))
rows = [rows[0]] + [[r_ if r_ else "" for r_ in row] for row in rows[1:]]
md_table_to_docx(doc, rows, font=8)
out4 = os.path.join(SUP, "STROBE_Checklist.docx")
doc.save(out4)
print("saved:", out3)
print("saved:", out4)
