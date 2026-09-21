"""Rebuild manuscript PDF bundles from deposited canonical individual pages."""
from pathlib import Path
from pypdf import PdfReader,PdfWriter

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/figure_bundles';OUT.mkdir(parents=True,exist_ok=True)
for stem,folder,pattern,n in [('Main_Figures_1_to_7','main','Figure_{:02d}.pdf',7),('Supplementary_Figures_1_to_12','supplementary','Supplementary_Figure_{:02d}.pdf',12),('Extended_Data_Figures_1_to_10','extended_data','Extended_Data_Figure_{:02d}.pdf',10)]:
    writer=PdfWriter()
    for i in range(1,n+1):
        source=ROOT/'figures'/folder/pattern.format(i)
        reader=PdfReader(source);assert len(reader.pages)==1,source
        writer.add_page(reader.pages[0])
    output=OUT/(stem+'.pdf');writer.write(output)
    assert len(PdfReader(output).pages)==n
    print(output)
