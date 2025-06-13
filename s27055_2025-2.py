import sys
from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt
Entrez.email = sys.argv[1]
Entrez.api_key = sys.argv[2]
taxid = sys.argv[3]
minL, maxL = map(int, sys.argv[4:6])
handle = Entrez.esearch(db="nucleotide", term=f"txid{taxid}[Orgn]", retmax=200)
ids = Entrez.read(handle)["IdList"]
stream = Entrez.efetch(db="nucleotide", id=",".join(ids),rettype="gb", retmode="text")
records = SeqIO.parse(stream, "genbank")
data = []
for rec in records:
    length = len(rec.seq)
    if minL <= length <= maxL:
        data.append((rec.id, length, rec.description))
df = pd.DataFrame(data, columns=["accession","length","description"])
df.to_csv("wyniki.csv", index=False)
df = df.sort_values("length", ascending=False)
plt.figure(figsize=(8,6))
plt.bar(df["accession"], df["length"])
plt.xticks(rotation=90)
plt.ylabel("Length")
plt.xlabel("Accession")
plt.tight_layout()
plt.savefig("wykres.png")