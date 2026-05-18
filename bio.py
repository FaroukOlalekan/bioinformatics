# ==========================================
# 🧬 SIMPLE INDUSTRY-STYLE BIOINFORMATICS PIPELINE
# DNA + RNA + PROTEIN ANALYSIS
# ==========================================

import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from collections import Counter
import math

# ---------------- PAGE ----------------
st.set_page_config(page_title="Bioinformatics Pipeline", layout="wide")

st.title("🧬 Bioinformatics Pipeline System")
st.write("Simple DNA, RNA, and Protein analysis pipeline")



# ==========================================
# CODON TABLE
# ==========================================

CODON_TABLE = {
    'ATA':'I','ATC':'I','ATT':'I','ATG':'M',
    'ACA':'T','ACC':'T','ACG':'T','ACT':'T',
    'AAC':'N','AAT':'N','AAA':'K','AAG':'K',
    'AGC':'S','AGT':'S','AGA':'R','AGG':'R',
    'CTA':'L','CTC':'L','CTG':'L','CTT':'L',
    'CCA':'P','CCC':'P','CCG':'P','CCT':'P',
    'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q',
    'CGA':'R','CGC':'R','CGG':'R','CGT':'R',
    'GTA':'V','GTC':'V','GTG':'V','GTT':'V',
    'GCA':'A','GCC':'A','GCG':'A','GCT':'A',
    'GAC':'D','GAT':'D','GAA':'E','GAG':'E',
    'GGA':'G','GGC':'G','GGG':'G','GGT':'G',
    'TCA':'S','TCC':'S','TCG':'S','TCT':'S',
    'TTC':'F','TTT':'F','TTA':'L','TTG':'L',
    'TAC':'Y','TAT':'Y','TAA':'*','TAG':'*',
    'TGC':'C','TGT':'C','TGA':'*','TGG':'W'
}

STOP_CODONS = ["TAA", "TAG", "TGA"]


# ==========================================
# FUNCTIONS
# ==========================================

def detect_type(seq):
    seq = seq.upper()

    if "U" in seq and "T" not in seq:
        return "RNA"

    elif set(seq).issubset(set("ATGCN")):
        return "DNA"

    return "PROTEIN"



def gc_content(seq):
    seq = seq.upper()

    if len(seq) == 0:
        return 0

    gc = seq.count("G") + seq.count("C")
    return round(gc / len(seq) * 100, 2)



def transcribe(seq):
    return seq.replace("T", "U")



def translate_dna(seq):
    protein = ""

    for i in range(0, len(seq)-2, 3):
        codon = seq[i:i+3]
        aa = CODON_TABLE.get(codon, "X")

        if aa == "*":
            break

        protein += aa

    return protein



def find_orfs(seq):
    orfs = []

    for frame in range(3):
        for i in range(frame, len(seq)-2, 3):
            codon = seq[i:i+3]

            if codon == "ATG":
                protein = translate_dna(seq[i:])

                if len(protein) > 10:
                    orfs.append(protein)

    return orfs



def entropy(seq):
    counts = Counter(seq)
    total = len(seq)

    if total == 0:
        return 0

    score = 0

    for base in counts:
        p = counts[base] / total

        if p > 0:
            score -= p * math.log2(p)

    return round(score, 2)


# ==========================================
# LOAD FASTA
# ==========================================

def load_fasta(file):
    ids = []
    sequences = []

    for record in SeqIO.parse(file, "fasta"):
        ids.append(record.id)
        sequences.append(str(record.seq).upper())

    return ids, sequences


# ==========================================
# PIPELINE ENGINE
# ==========================================

class Pipeline:

    def __init__(self, sequences):
        self.sequences = sequences

    def qc(self):
        return {
            "Total Sequences": len(self.sequences),
            "Average Length": round(sum(len(s) for s in self.sequences) / len(self.sequences), 2),
            "Max Length": max(len(s) for s in self.sequences),
            "Min Length": min(len(s) for s in self.sequences)
        }

    def gc_analysis(self):
        values = [gc_content(s) for s in self.sequences if detect_type(s) != "PROTEIN"]

        return {
            "Average GC": round(sum(values)/len(values), 2) if values else 0,
            "Values": values
        }

    def orf_analysis(self):
        results = []

        for seq in self.sequences:
            if detect_type(seq) == "DNA":
                results.extend(find_orfs(seq))

        return results


# ==========================================
# UI
# ==========================================

uploaded = st.file_uploader("Upload FASTA file", type=["fasta", "fa", "txt"])

question = st.text_input(
    "Ask pipeline question",
    placeholder="GC? ORF? longest? protein?"
)

if uploaded:

    ids, sequences = load_fasta(uploaded)

    pipeline = Pipeline(sequences)

    # ---------------- TABLE ----------------
    st.subheader("📄 Sequence Summary")

    table = pd.DataFrame({
        "ID": ids,
        "Type": [detect_type(s) for s in sequences],
        "Length": [len(s) for s in sequences],
        "GC%": [gc_content(s) if detect_type(s) != "PROTEIN" else None for s in sequences]
    })

    st.dataframe(table)

    # ---------------- QC ----------------
    st.subheader("🧪 QC Metrics")
    qc = pipeline.qc()
    st.json(qc)

    # ---------------- GC ----------------
    st.subheader("🧬 GC Analysis")
    gc = pipeline.gc_analysis()
    st.json(gc)

    fig, ax = plt.subplots()
    ax.hist(gc["Values"], bins=10)
    ax.set_xlabel("GC %")
    ax.set_ylabel("Frequency")
    st.pyplot(fig)

    # ---------------- ORF ----------------
    st.subheader("🧬 ORF Detection")
    orfs = pipeline.orf_analysis()

    if orfs:
        st.write(f"Detected {len(orfs)} ORFs")
        st.code(orfs[:5])
    else:
        st.write("No ORFs detected")

    # ---------------- Q&A ----------------
    if question:

        q = question.lower()

        st.subheader("🤖 Pipeline Answer")

        if "gc" in q:
            st.write(gc)

        elif "orf" in q:
            st.write(f"Detected {len(orfs)} ORFs")

        elif "protein" in q:
            proteins = [translate_dna(s) for s in sequences if detect_type(s) == "DNA"]
            st.write(proteins[:5])

        elif "longest" in q:
            longest = max(sequences, key=len)
            st.write(f"Longest sequence length = {len(longest)}")

        else:
            st.write("Question not recognized")

else:
    st.info("Upload FASTA sequences to start pipeline analysis")


# ==========================================
# IMPORTANT NOTE
# ==========================================

st.markdown("""
### ⚠️ Simplified Pipeline

This version now includes:

- DNA analysis
- RNA analysis
- Protein translation
- ORF detection
- GC content
- QC metrics
- Simple pipeline Q&A

But it is still simplified:

- no BLAST integration
- no FASTQ quality scores
- no alignment tools
- no Snakemake workflow
- no real HPC execution

This is a strong learning and portfolio architecture.
""")
