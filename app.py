import streamlit as st
import pandas as pd

from database import load_fasta_database
from search_engine import similarity_search


# =========================================================
# PAGE CONFIG
# =========================================================

st.set_page_config(
    page_title="BioBLAST-X",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)


# =========================================================
# CUSTOM CSS
# =========================================================

st.markdown("""
<style>

.main {
    background-color: #0e1117;
}

.block-container {
    padding-top: 2rem;
    padding-bottom: 3rem;
}

.hero {
    padding: 25px;
    border-radius: 15px;
    background: linear-gradient(
        135deg,
        #172554,
        #0f766e
    );
    margin-bottom: 25px;
}

.hero h1 {
    color: white;
    font-size: 42px;
    margin-bottom: 5px;
}

.hero p {
    color: #dbeafe;
    font-size: 17px;
}

.section-title {
    font-size: 26px;
    font-weight: 700;
    margin-top: 25px;
    margin-bottom: 15px;
}

.metric-card {
    padding: 18px;
    border-radius: 12px;
    background-color: #171b26;
    border: 1px solid #2b3242;
    text-align: center;
}

.metric-value {
    font-size: 28px;
    font-weight: bold;
}

.metric-label {
    color: #9ca3af;
    font-size: 14px;
}

.result-card {
    padding: 20px;
    border-radius: 14px;
    background-color: #171b26;
    border: 1px solid #2b3242;
    margin-bottom: 15px;
}

</style>
""", unsafe_allow_html=True)


# =========================================================
# HEADER
# =========================================================

st.markdown("""
<div class="hero">

<h1>🧬 BioBLAST-X</h1>

<p>
BLAST-Inspired Biological Sequence Similarity Analyzer
</p>

<p>
Analyze DNA sequences using k-mer based similarity screening
and sequence comparison.
</p>

</div>
""", unsafe_allow_html=True)


# =========================================================
# SIDEBAR
# =========================================================

st.sidebar.title("⚙️ Analysis Settings")

st.sidebar.markdown("### Search Parameters")

k = st.sidebar.slider(
    "K-mer Size",
    min_value=2,
    max_value=8,
    value=4
)

st.sidebar.markdown("---")

st.sidebar.markdown("""
### 🧬 Pipeline

**1. Sequence Input**

↓

**2. Sequence Validation**

↓

**3. K-mer Screening**

↓

**4. Similarity Search**

↓

**5. Identity Calculation**

↓

**6. Similarity Ranking**
""")


# =========================================================
# QUERY INPUT
# =========================================================

st.markdown(
    '<div class="section-title">1️⃣ Query Sequence</div>',
    unsafe_allow_html=True
)

input_method = st.radio(
    "Choose input method",
    ["Paste Sequence", "Upload FASTA File"],
    horizontal=True
)

query_sequence = ""


# =========================================================
# PASTE SEQUENCE
# =========================================================

if input_method == "Paste Sequence":

    query_sequence = st.text_area(
        "Enter DNA sequence",
        height=160,
        placeholder=(
            "Example:\n"
            "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
        )
    )


# =========================================================
# FASTA UPLOAD
# =========================================================

else:

    uploaded_file = st.file_uploader(
        "Upload FASTA file",
        type=["fasta", "fa", "txt"]
    )

    if uploaded_file is not None:

        content = uploaded_file.read().decode(
            "utf-8"
        )

        sequence_lines = []

        for line in content.splitlines():

            line = line.strip()

            if line and not line.startswith(">"):
                sequence_lines.append(line)

        query_sequence = "".join(
            sequence_lines
        )

        st.success(
            f"✅ FASTA file loaded: {uploaded_file.name}"
        )


# =========================================================
# CLEAN SEQUENCE
# =========================================================

query_sequence = (
    query_sequence
    .upper()
    .replace(" ", "")
    .replace("\n", "")
    .replace("\r", "")
    .replace("\t", "")
)


# =========================================================
# VALIDATION
# =========================================================

if query_sequence:

    valid_bases = set("ATGCN")

    invalid_characters = sorted(
        set(query_sequence) - valid_bases
    )

    if invalid_characters:

        st.error(
            "❌ Invalid DNA sequence detected."
        )

        st.write(
            "Invalid characters:"
        )

        st.code(
            ", ".join(invalid_characters)
        )

        st.stop()

    else:

        st.success(
            "✅ Sequence is valid."
        )


        # =================================================
        # SEQUENCE STATISTICS
        # =================================================

        st.markdown(
            '<div class="section-title">2️⃣ Sequence Statistics</div>',
            unsafe_allow_html=True
        )

        sequence_length = len(query_sequence)

        gc_count = (
            query_sequence.count("G")
            +
            query_sequence.count("C")
        )

        gc_percentage = (
            gc_count / sequence_length
        ) * 100


        col1, col2, col3 = st.columns(3)


        with col1:

            st.markdown(
                f"""
                <div class="metric-card">

                <div class="metric-value">
                {sequence_length}
                </div>

                <div class="metric-label">
                Sequence Length (bp)
                </div>

                </div>
                """,
                unsafe_allow_html=True
            )


        with col2:

            st.markdown(
                f"""
                <div class="metric-card">

                <div class="metric-value">
                {gc_percentage:.2f}%
                </div>

                <div class="metric-label">
                GC Content
                </div>

                </div>
                """,
                unsafe_allow_html=True
            )


        with col3:

            st.markdown(
                f"""
                <div class="metric-card">

                <div class="metric-value">
                {k}
                </div>

                <div class="metric-label">
                K-mer Size
                </div>

                </div>
                """,
                unsafe_allow_html=True
            )


        # =================================================
        # PROCESSED SEQUENCE
        # =================================================

        st.markdown(
            '<div class="section-title">3️⃣ Processed Sequence</div>',
            unsafe_allow_html=True
        )

        st.code(
            query_sequence,
            language="text"
        )


        # =================================================
        # RUN ANALYSIS
        # =================================================

        st.markdown(
            '<div class="section-title">4️⃣ Similarity Analysis</div>',
            unsafe_allow_html=True
        )

        run_analysis = st.button(
            "🔍 Run BioBLAST-X",
            type="primary",
            use_container_width=True
        )


        if run_analysis:

            try:

                with st.spinner(
                    "Running sequence similarity analysis..."
                ):

                    # Load local FASTA database

                    database = load_fasta_database(
                        "data/sequences.fasta"
                    )

                    # Run similarity search

                    results = similarity_search(
                        query_sequence,
                        database,
                        k=k
                    )


                # =================================================
                # SUCCESS
                # =================================================

                st.success(
                    f"✅ Analysis completed successfully. "
                    f"{len(results)} sequences analyzed."
                )


                # =================================================
                # RESULTS
                # =================================================

                st.markdown(
                    '<div class="section-title">5️⃣ Similarity Results</div>',
                    unsafe_allow_html=True
                )


                result_rows = []


                for rank, result in enumerate(
                    results,
                    start=1
                ):

                    result_rows.append({

                        "Rank": rank,

                        "Sequence ID":
                            result.get("id", "Unknown"),

                        "Identity (%)":
                            result.get("identity", 0),

                        "Shared K-mers":
                            result.get("shared_kmers", 0)

                    })


                results_df = pd.DataFrame(
                    result_rows
                )


                # =================================================
                # RESULTS TABLE
                # =================================================

                st.dataframe(
                    results_df,
                    use_container_width=True,
                    hide_index=True
                )


                # =================================================
                # SIMILARITY GRAPH
                # =================================================

                st.markdown(
                    '<div class="section-title">6️⃣ Similarity Visualization</div>',
                    unsafe_allow_html=True
                )


                chart_df = results_df[
                    ["Sequence ID", "Identity (%)"]
                ].copy()

                chart_df = chart_df.set_index(
                    "Sequence ID"
                )

                st.bar_chart(
                    chart_df,
                    y="Identity (%)"
                )


                # =================================================
                # BEST MATCH
                # =================================================

                if len(results) > 0:

                    best_hit = results[0]


                    st.markdown(
                        '<div class="section-title">7️⃣ Top Similarity Match</div>',
                        unsafe_allow_html=True
                    )


                    best_id = best_hit.get(
                        "id",
                        "Unknown"
                    )

                    best_identity = best_hit.get(
                        "identity",
                        0
                    )

                    best_kmers = best_hit.get(
                        "shared_kmers",
                        0
                    )


                    col1, col2, col3 = st.columns(3)


                    with col1:

                        st.metric(
                            "Top Match",
                            best_id
                        )


                    with col2:

                        st.metric(
                            "Identity",
                            f"{best_identity:.2f}%"
                        )


                    with col3:

                        st.metric(
                            "Shared K-mers",
                            best_kmers
                        )


                    # =============================================
                    # OPTIONAL ALIGNMENT INFORMATION
                    # =============================================

                    if (
                        "aligned_query"
                        in best_hit
                        and
                        "aligned_target"
                        in best_hit
                    ):

                        st.markdown(
                            '<div class="section-title">8️⃣ Sequence Alignment</div>',
                            unsafe_allow_html=True
                        )

                        st.write(
                            f"**Query:** {best_id}"
                        )

                        st.code(
                            f"Query  : "
                            f"{best_hit['aligned_query']}\n\n"
                            f"Target : "
                            f"{best_hit['aligned_target']}"
                        )


                # =================================================
                # DOWNLOAD RESULTS
                # =================================================

                st.markdown(
                    '<div class="section-title">9️⃣ Export Results</div>',
                    unsafe_allow_html=True
                )


                csv_data = results_df.to_csv(
                    index=False
                )


                st.download_button(
                    label="⬇️ Download Results as CSV",
                    data=csv_data,
                    file_name="bioblast_x_results.csv",
                    mime="text/csv",
                    use_container_width=True
                )


            except Exception as e:

                st.error(
                    "❌ An error occurred during analysis."
                )

                st.code(
                    str(e)
                )


# =========================================================
# EMPTY STATE
# =========================================================

else:

    st.info(
        "🧬 Enter a DNA sequence or upload a FASTA file "
        "to begin analysis."
    )


# =========================================================
# FOOTER
# =========================================================

st.markdown("---")

st.caption(
    "BioBLAST-X | Python • Streamlit • Bioinformatics"
)