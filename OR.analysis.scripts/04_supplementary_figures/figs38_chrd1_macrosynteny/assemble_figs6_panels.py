import argparse
from pathlib import Path

from PIL import Image
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf._page import PageObject


parser = argparse.ArgumentParser(description="Assemble Supplementary Figure S6 panels")
parser.add_argument("--panel-dir", type=Path, default=Path(__file__).resolve().parent)
ROOT = parser.parse_args().panel_dir.resolve()
PDF_OUT = ROOT / "Supplementary_Figure_S6_combined.pdf"
PNG_OUT = ROOT / "Supplementary_Figure_S6_combined.png"
PDF_A = ROOT / "panel_a_macro_alignment.pdf"
PDF_B = ROOT / "panel_b_regional_homology.pdf"
PNG_A = ROOT / "panel_a_macro_alignment.png"
PNG_B = ROOT / "panel_b_regional_homology.png"

# A modest white divider: 2 mm in PDF points and 24 pixels at 300 dpi.
GAP_PT = 72 / 25.4 * 2
GAP_PX = 24


def combine_pdf() -> tuple[float, float]:
    a = PdfReader(str(PDF_A)).pages[0]
    b = PdfReader(str(PDF_B)).pages[0]
    aw, ah = float(a.mediabox.width), float(a.mediabox.height)
    bw, bh = float(b.mediabox.width), float(b.mediabox.height)
    if abs(aw - bw) > 0.1:
        raise ValueError(f"Panel widths differ: A={aw:.2f} pt, B={bw:.2f} pt")

    page = PageObject.create_blank_page(width=aw, height=ah + GAP_PT + bh)
    page.merge_transformed_page(a, Transformation().translate(tx=0, ty=bh + GAP_PT))
    page.merge_transformed_page(b, Transformation().translate(tx=0, ty=0))
    writer = PdfWriter()
    writer.add_page(page)
    with PDF_OUT.open("wb") as handle:
        writer.write(handle)
    return aw, ah + GAP_PT + bh


def combine_png() -> tuple[int, int]:
    with Image.open(PNG_A) as top, Image.open(PNG_B) as bottom:
        if top.width != bottom.width:
            raise ValueError(f"Panel widths differ: A={top.width}px, B={bottom.width}px")
        canvas = Image.new("RGB", (top.width, top.height + GAP_PX + bottom.height), "white")
        canvas.paste(top.convert("RGB"), (0, 0))
        canvas.paste(bottom.convert("RGB"), (0, top.height + GAP_PX))
        canvas.save(PNG_OUT, dpi=(300, 300), optimize=True)
        return canvas.size


if __name__ == "__main__":
    pdf_w, pdf_h = combine_pdf()
    png_w, png_h = combine_png()
    print(f"PDF: {PDF_OUT} ({pdf_w:.2f} x {pdf_h:.2f} pt)")
    print(f"PNG: {PNG_OUT} ({png_w} x {png_h} px, 300 dpi)")
