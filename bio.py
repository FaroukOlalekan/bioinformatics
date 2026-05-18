import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from collections import Counter
from io import StringIO
import math

# ======================================================
# PAGE CONFIG
# ======================================================

st.set_page_config(
    page_title="Bioinformatics Pipeline System",
    layout="wide"
)

# ======================================================
# CUSTOM STYLING
# ======================================================

st.markdown("""
<style>

.stApp {
    background: linear-gradient(to right, #0f172a, #1e293b);
    color: #f8fafc;
}

h1, h2, h3, h4 {
    color: #38bdf8;
}

p, label, div {
    color: #e2e8f0;
}

[data-testid="stSidebar"] {
    background-color: #111827;
}

.stButton>button {
    background-color: #38bdf8;
    color: black;
    border-radius: 10px;
    border: none;
    padding: 0.5rem 1rem;
    font-weight: bold;
}

.stButton>button:hover {
    background-color: #0ea5e9;
    color: white;
}

.stTextInput>div>div>input {
    background-color: #1e293b;
    color: white;
}

.stFileUploader {
    background-color: #1e293b;
    border-radius: 10px;
    padding: 10px;
}

</style>
""", unsafe_allow_html=True)

# ======================================================
# TITLE
# ======================================================

st.title("🧬 Multi-Omics Bioinformatics Pipeline")

st.write("""
Analyze DNA, RNA, and Protein FASTA files with:

- QC Metrics
- GC Content Analysis
- ORF Detection
- DNA Translation
- Entropy Analysis
- BLAST Search
- Interactive Bioinformatics Questions
""")



# ======================================================
# CODON TABLE
# ======================================================

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

# ======================================================
# FUNCTIONS
# ======================================================

def detect_type(seq):

    seq = seq.upper()

    if "U" in seq and "T" not in seq:
        return "RNA"

    elif set(seq).issubset(set("ATGCN")):
        return "DNA"

    else:
        return "PROTEIN"


def gc_content(seq):

    seq = seq.upper()

    if len(seq) == 0:
        return 0

    gc = seq.count("G") + seq.count("C")

    return round((gc / len(seq)) * 100, 2)


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


# ======================================================
# FASTA LOADER
# ======================================================

def load_fasta(file):

    ids = []
    sequences = []

    try:

        fasta_text = StringIO(
            file.getvalue().decode("utf-8")
        )

        for record in SeqIO.parse(fasta_text, "fasta"):

            ids.append(record.id)
            sequences.append(str(record.seq).upper())

        return ids, sequences

    except Exception as e:

        st.error(f"FASTA parsing error: {e}")

        return [], []


# ======================================================
# BLAST FUNCTION
# ======================================================

def run_blast(sequence_id, sequence):

    try:

        sequence_type = detect_type(sequence)

        # DNA
        if sequence_type == "DNA":

            program = "blastn"
            database = "nt"

        # RNA
        elif sequence_type == "RNA":

            program = "blastn"
            database = "nt"

            sequence = sequence.replace("U", "T")

        # PROTEIN
        else:

            program = "blastp"
            database = "nr"

        result_handle = NCBIWWW.qblast(
            program,
            database,
            sequence
        )

        blast_record = NCBIXML.read(result_handle)

        results = []

        for alignment in blast_record.alignments[:10]:

            for hsp in alignment.hsps[:1]:

                results.append({

                    "Query": sequence_id,
                    "Type": sequence_type,
                    "Title": alignment.title,
                    "Length": alignment.length,
                    "E-Value": hsp.expect,
                    "Score": hsp.score

                })

        return results

    except Exception as e:

        st.error(f"BLAST Error: {e}")

        return []


# ======================================================
# PIPELINE CLASS
# ======================================================

class Pipeline:

    def __init__(self, sequences):

        self.sequences = sequences

    def qc(self):

        return {

            "Total Sequences": len(self.sequences),

            "Average Length": round(
                sum(len(s) for s in self.sequences) / len(self.sequences),
                2
            ),

            "Max Length": max(len(s) for s in self.sequences),

            "Min Length": min(len(s) for s in self.sequences)

        }

    def gc_analysis(self):

        values = [

            gc_content(s)

            for s in self.sequences

            if detect_type(s) != "PROTEIN"

        ]

        return {

            "Average GC": round(
                sum(values) / len(values),
                2
            ) if values else 0,

            "Values": values

        }

    def orf_analysis(self):

        results = []

        for seq in self.sequences:

            if detect_type(seq) == "DNA":

                results.extend(find_orfs(seq))

        return results


# ======================================================
# FILE UPLOAD
# ======================================================

uploaded = st.file_uploader(
    "📂 Upload FASTA File",
    type=["fasta", "fa", "faa", "fna"]
)

question = st.text_input(
    "❓ Ask Pipeline Question",
    placeholder="GC analysis? ORFs? Longest sequence?"
)

# ======================================================
# MAIN APPLICATION
# ======================================================

if uploaded:

    ids, sequences = load_fasta(uploaded)

    if sequences:

        pipeline = Pipeline(sequences)

        # ==================================================
        # SUMMARY TABLE
        # ==================================================

        st.subheader("📄 Sequence Summary")

        table = pd.DataFrame({

            "ID": ids,

            "Type": [
                detect_type(s)
                for s in sequences
            ],

            "Length": [
                len(s)
                for s in sequences
            ],

            "GC%": [
                gc_content(s)
                if detect_type(s) != "PROTEIN"
                else None
                for s in sequences
            ],

            "Entropy": [
                entropy(s)
                for s in sequences
            ]

        })

        st.dataframe(table)

        # ==================================================
        # QC
        # ==================================================

        st.subheader("🧪 QC Metrics")

        qc = pipeline.qc()

        st.json(qc)

        # ==================================================
        # GC ANALYSIS
        # ==================================================

        st.subheader("🧬 GC Content Analysis")

        gc = pipeline.gc_analysis()

        st.json(gc)

        fig, ax = plt.subplots()

        ax.hist(gc["Values"], bins=10)

        ax.set_xlabel("GC %")
        ax.set_ylabel("Frequency")

        st.pyplot(fig)

        # ==================================================
        # ORF DETECTION
        # ==================================================

        st.subheader("🧬 ORF Detection")

        orfs = pipeline.orf_analysis()

        if orfs:

            st.success(f"Detected {len(orfs)} ORFs")

            st.code(orfs[:5])

        else:

            st.warning("No ORFs detected")

        # ==================================================
        # PROTEIN TRANSLATION
        # ==================================================

        st.subheader("🧬 Protein Translation")

        proteins = [

            translate_dna(seq)

            for seq in sequences

            if detect_type(seq) == "DNA"

        ]

        if proteins:

            st.code(proteins[:5])

        # ==================================================
        # BLAST
        # ==================================================

        st.subheader("🧬 Multi-Omics BLAST Search")

        blast_option = st.checkbox(
            "Run BLAST for DNA, RNA, and Protein sequences"
        )

        if blast_option:

            all_results = []

            sequence_data = [

                (ids[i], sequences[i])

                for i in range(len(sequences))

            ]

            with st.spinner(
                "Running BLAST for all sequences..."
            ):

                for seq_id, seq in sequence_data:

                    results = run_blast(
                        seq_id,
                        seq
                    )

                    all_results.extend(results)

            if all_results:

                blast_df = pd.DataFrame(all_results)

                st.dataframe(blast_df)

            else:

                st.warning(
                    "No BLAST results found"
                )

        # ==================================================
        # AI-LIKE QUESTIONS
        # ==================================================

        if question:

            q = question.lower()

            st.subheader("🤖 Pipeline Answer")

            if "gc" in q:

                st.write(gc)

            elif "orf" in q:

                st.write(f"Detected {len(orfs)} ORFs")

            elif "protein" in q:

                st.write(proteins[:5])

            elif "longest" in q:

                longest = max(sequences, key=len)

                st.write(
                    f"Longest sequence length = {len(longest)}"
                )

            elif "entropy" in q:

                st.write(
                    table[["ID", "Entropy"]]
                )

            else:

                st.write(
                    "Question not recognized."
                )

else:

    st.info(
        "Upload FASTA sequences to begin analysis."
    )

# ======================================================
# FOOTER
# ======================================================

st.markdown("""
---
### 👨‍💻 Developed by Farouk

Multi-Omics Bioinformatics Pipeline System
""")
