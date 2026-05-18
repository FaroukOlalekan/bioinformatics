# 🧬 Bioinformatics Pipeline System

A simple multi-omics bioinformatics pipeline application built with Python and Streamlit.

This system supports:

- DNA analysis
- RNA analysis
- Protein analysis
- GC content analysis
- ORF detection
- Sequence quality control (QC)
- Protein translation
- FASTA file parsing
- Interactive bioinformatics Q&A

---

# 🚀 Features

## ✅ DNA Analysis
- GC content calculation
- ORF detection
- DNA → RNA transcription
- DNA → Protein translation

## ✅ RNA Analysis
- RNA sequence identification
- RNA translation support

## ✅ Protein Analysis
- Protein sequence recognition
- Amino acid analysis support

## ✅ QC Metrics
- Sequence length statistics
- Maximum sequence length
- Minimum sequence length
- Average sequence length

## ✅ Visualization
- GC distribution plots
- Interactive tables

## ✅ Interactive Pipeline Questions
Users can ask questions like:

- "Show GC analysis"
- "Find ORFs"
- "Show protein sequences"
- "Longest sequence?"

---

# 📂 Supported File Types

- `.fasta`
- `.fa`
- `.txt`

---

# 🛠️ Installation

## 1. Clone Repository

```bash
git clone https://github.com/YOUR_USERNAME/bioinformatics-pipeline-system.git
cd bioinformatics-pipeline-system
```

---

## 2. Install Dependencies

```bash
pip install -r requirements.txt
```

---

# ▶️ Run Application

```bash
streamlit run app.py
```

The application will open in your browser:

```text
http://localhost:8501
```

---

# 📊 Technologies Used

- Python
- Streamlit
- Pandas
- Matplotlib
- Biopython

---

# 🧠 System Architecture

```text
User Upload
     ↓
FASTA Parser
     ↓
Pipeline Engine
     ↓
DNA / RNA / Protein Analysis
     ↓
QC + ORF + Translation
     ↓
Visualization + Q&A
```

---

# ⚠️ Current Limitations

This is a simplified educational bioinformatics pipeline.

Current limitations:

- No BLAST integration
- No FASTQ quality scores
- No sequence alignment tools
- No Snakemake workflow
- No HPC support
- No database integration

---

# 🚀 Future Improvements

Planned upgrades:

- BLAST integration
- FastQC support
- MultiQC reporting
- Sequence alignment
- Variant calling
- Docker support
- Snakemake workflows
- AI-powered biological interpretation

---

# 👨‍💻 Author

Farouk

MSc Bioinformatics Student

---

# 📜 License

MIT License