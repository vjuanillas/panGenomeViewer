import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from io import StringIO

# Page configuration
st.set_page_config(
    page_title="Pangenome Statistics Explorer",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom UI Styling
st.markdown("""
    <style>
    .main { background-color: #f9fbfc; }
    .stMetric {
        background-color: #ffffff;
        padding: 20px;
        border-radius: 12px;
        box-shadow: 0 4px 6px rgba(0,0,0,0.03);
        border-top: 4px solid #3b82f6;
    }
    .compartment-card {
        padding: 20px;
        border-radius: 10px;
        color: white;
        text-align: center;
        margin-bottom: 15px;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
    }
    .stTabs [data-baseweb="tab"] {
        font-weight: 600;
        padding: 10px 20px;
    }
    </style>
    """, unsafe_allow_html=True)

# --- CACHED DATA PARSING ---

@st.cache_data(show_spinner=False)
def process_file(file_contents, file_name):
    try:
        if len(file_contents) == 0:
            return None, "Empty file"

        lines = file_contents.decode("utf-8").splitlines()
        is_panacus = any("panacus" in line.lower() for line in lines[:5])
        
        data_lines = [l for l in lines if not l.startswith("#")]
        if not data_lines:
            return None, "No data rows found"

        if is_panacus:
            # 1. Panacus Growth check (e.g. ordered-histgrowth)
            coverage_line = next((l for l in lines if l.startswith("coverage")), None)
            if coverage_line or any("ordered-growth" in l for l in lines[:10]):
                # Use the 'coverage' line to get column headers (1, 2, 5, 15 etc)
                if coverage_line:
                    cols = ['Sample'] + coverage_line.split('\t')[1:]
                else:
                    cols = None
                
                # Find where actual data rows begin
                start_idx = 0
                for i, l in enumerate(data_lines):
                    parts = l.split('\t')
                    # Data rows usually start with a string (sample name) then a number
                    if len(parts) > 1 and parts[1].strip().replace('.','',1).isdigit():
                        start_idx = i
                        break
                
                df = pd.read_csv(StringIO("\n".join(data_lines[start_idx:])), sep='\t', names=cols)
                return df, "panacus_growth"

            # 2. Panacus Info / Table style (feature category countable value)
            first_row_split = data_lines[0].split('\t')
            if "feature" in first_row_split and "category" in first_row_split:
                df = pd.read_csv(StringIO("\n".join(data_lines)), sep='\t')
                return df, "panacus_table"

        # 3. Generic/PA Matrix fallback
        df = pd.read_csv(StringIO("\n".join(data_lines)), sep='\t')
        if len(df.columns) > 2:
            numeric_cols = df.select_dtypes(include=[np.number]).columns
            if len(numeric_cols) >= (len(df.columns) - 1):
                return df, "pa_matrix"

        return df, "generic_stats"
    except Exception as e:
        return None, str(e)

# --- SIDEBAR ---

with st.sidebar:
    st.title("Pangenome Statistics Explorer")
    st.divider()
    
    st.header("Data Upload")
    uploaded_files = st.file_uploader(
        "Upload Stats/Tables/Matrices", 
        type=["tsv", "csv", "txt"], 
        accept_multiple_files=True
    )

    if 'detected_samples' not in st.session_state:
        st.session_state.detected_samples = []

    st.divider()
    st.header("🧪 Sample Filtering")
    exclude_samples = st.multiselect(
        "Exclude Samples",
        options=st.session_state.detected_samples,
        help="Exclude outgroups. Matrix and growth calculations update instantly."
    )

# --- PROCESS FILES ---

registry = {}
new_samples = set()

if uploaded_files:
    with st.spinner("Processing files..."):
        for f in uploaded_files:
            file_bytes = f.getvalue()
            result, dtype_or_error = process_file(file_bytes, f.name)
            if result is not None:
                registry[f.name] = {"df": result, "type": dtype_or_error}
                
                if dtype_or_error == "panacus_growth":
                    new_samples.update(result['Sample'].unique().tolist())
                elif dtype_or_error == "panacus_table":
                    # Samples are found in 'category' column when 'feature' is 'group'
                    s_list = result[result['feature'] == 'group']['category'].unique().tolist()
                    new_samples.update(s_list)
                elif dtype_or_error == "pa_matrix":
                    new_samples.update(result.columns[1:].tolist())
            else:
                st.error(f"Error in {f.name}: {dtype_or_error}")

# Sort and sync detected samples
sorted_samples = sorted(list(new_samples))
if sorted_samples != st.session_state.detected_samples:
    st.session_state.detected_samples = sorted_samples
    st.rerun()

# --- MAIN DASHBOARD ---

if not registry:
    st.title("Pangenome Analytics")
    st.info("Upload Panacus growth (.tsv) or info (.tsv) files to begin.")
else:
    active_registry = {}
    for name, data in registry.items():
        df = data['df'].copy()
        if data['type'] == "panacus_growth":
            df = df[~df['Sample'].isin(exclude_samples)].reset_index(drop=True)
        elif data['type'] == "panacus_table":
            # Filter rows where category is an excluded sample
            mask = ~((df['feature'] == 'group') & (df['category'].isin(exclude_samples)))
            df = df[mask]
        elif data['type'] == "pa_matrix":
            cols_to_drop = [c for c in exclude_samples if c in df.columns]
            df = df.drop(columns=cols_to_drop)
        active_registry[name] = {"df": df, "type": data['type']}

    st.title("Pangenome Analysis Results Viewer")
    
    growth_files = [n for n, d in active_registry.items() if d['type'] == "panacus_growth"]
    table_files = [n for n, d in active_registry.items() if d['type'] == "panacus_table"]
    matrix_files = [n for n, d in active_registry.items() if d['type'] == "pa_matrix"]

    t_overview, t_growth, t_nodes = st.tabs(["Graph Overview", "Ordered Growth", "Node Categories"])

    with t_overview:
        if table_files:
            df_t = active_registry[table_files[0]]['df']
            totals = df_t[df_t['category'] == 'total']
            if not totals.empty:
                m1, m2, m3, m4 = st.columns(4)
                
                def get_v(c):
                    r = totals[totals['countable'] == c]['value']
                    return int(r.iloc[0]) if not r.empty else 0
                
                m1.metric("Total Nodes", f"{get_v('node'):,}")
                m2.metric("Total BP", f"{get_v('bp'):,}")
                m3.metric("Total Edges", f"{get_v('edge'):,}")
                
                active_count = len(st.session_state.detected_samples) - len(exclude_samples)
                m4.metric("Active Samples", f"{active_count}")

            groups = df_t[df_t['feature'] == 'group']
            if not groups.empty:
                pivot = groups.pivot(index='category', columns='countable', values='value').reset_index()
                fig = px.bar(pivot, x='category', y='bp', title="Total length (basepairs) per Active Sample", color='bp')
                # Rename the x-axis and the y-axis labels
                fig.update_layout(
                    xaxis_title='Genomes',
                    yaxis_title='Total Length (bp)'
                )
                st.plotly_chart(fig, use_container_width=True)

    with t_growth:
        if growth_files:
            df_g = active_registry[growth_files[0]]['df']
            if not df_g.empty:
                # Filter out 'Sample' and any helper cols
                plot_cols = [c for c in df_g.columns if c not in ['Sample', 'Rank', 'Order']]
                df_g['Order'] = range(1, len(df_g) + 1)
                
                fig = go.Figure()
                for c in plot_cols:
                    fig.add_trace(go.Scatter(x=df_g['Order'], y=df_g[c], name=f"Shared by ≥ {c}", mode='lines+markers'))
                
                fig.update_layout(xaxis_title="Genomes Added", yaxis_title="Cumulative Length(basepair)", hovermode="x unified")
                st.plotly_chart(fig, use_container_width=True)

    with t_nodes:
        # Check if we have a growth file to calculate compartments (Core, Dispensable, Private)
        if growth_files:
            st.header("Pangenome Countable categories")
            df_c = active_registry[growth_files[0]]['df']
            
            if not df_c.empty:
                last_pt = df_c.iloc[-1]
                # Columns represent coverage thresholds
                cov_cols = [c for c in df_c.columns if c not in ['Sample', 'Order', 'Rank']]
                
                # total_bp is the 1x level
                total_pan = last_pt[cov_cols[0]]
                # shared_2 is the sequence in at least 2 genomes
                shared_2plus = last_pt[cov_cols[1]] if len(cov_cols) > 1 else total_pan
                
                # Core is sequence in ALL active genomes
                num_active = len(df_c)
                # Look for column matching current sample count, or highest available
                core_col = str(num_active) if str(num_active) in cov_cols else cov_cols[-1]
                core_seq = last_pt[core_col]
                
                private_seq = total_pan - shared_2plus
                dispensable_seq = shared_2plus - core_seq
                
                c1, c2, c3 = st.columns(3)
                c1.markdown(f'<div class="compartment-card" style="background-color: #166534;"><h4>Core</h4><h2>{core_seq/1e6:.2f} Mb</h2><p>Shared by {core_col} Genomes</p></div>', unsafe_allow_html=True)
                c2.markdown(f'<div class="compartment-card" style="background-color: #854d0e;"><h4>Dispensable</h4><h2>{dispensable_seq/1e6:.2f} Mb</h2><p>Shared by 2 to {int(core_col)-1}</p></div>', unsafe_allow_html=True)
                c3.markdown(f'<div class="compartment-card" style="background-color: #991b1b;"><h4>Private</h4><h2>{private_seq/1e6:.2f} Mb</h2><p>Unique to One</p></div>', unsafe_allow_html=True)
                
                fig_p = px.pie(
                    names=['Core', 'Dispensable', 'Private'],
                    values=[core_seq, dispensable_seq, private_seq],
                    color=['Core', 'Dispensable', 'Private'],
                    color_discrete_map={'Core':'#166534', 'Dispensable':'#854d0e', 'Private':'#991b1b'},
                    hole=0.4, title="Sequence Allocation (%)"
                )
                st.plotly_chart(fig_p, use_container_width=True)
            else:
                st.warning("No data remains after filtering.")
        else:
            st.warning("Upload an 'ordered-histgrowth' file to enable Compartment analysis.")

st.divider()
st.caption("Pangenome Explorer | Optimized for Panacus landraces results")