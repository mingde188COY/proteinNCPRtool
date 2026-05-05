import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ==========================================================
# Sidebar (LEFT instructions)
# ==========================================================
st.sidebar.title("Instructions")

st.sidebar.markdown("""
### Protein NCPR Heatmap Tool

This tool calculates the Net Charge Per Residue (NCPR) using a sliding window approach.

#### How to use:
1. Enter protein name
2. Paste amino acid sequence
3. Select sliding window size
4. Adjust visualization settings
5. Click **Generate**

#### Parameters:
- **Sliding window**: smoothing window for charge averaging
- **Tick interval**: spacing of residue labels on x-axis
- **Heatmap height**: vertical scaling of the plot
- **Font settings**: control text appearance

#### Output:
- NCPR heatmap (blue = negative, red = positive)
- Color scale bar
- Downloadable PNG figure
""")

# ==========================================================
# charge function
# ==========================================================
def get_charge(aa):
    if aa in ['K', 'R']:
        return 1
    elif aa in ['D', 'E']:
        return -1
    else:
        return 0

# ==========================================================
# NCPR calculation
# ==========================================================
def compute_ncpr(seq, window):
    charges = np.array([get_charge(a) for a in seq])
    half = window // 2

    ncpr = []
    for i in range(len(seq)):
        start = max(0, i - half)
        end = min(len(seq), i + half + 1)
        ncpr.append(np.mean(charges[start:end]))

    return np.array(ncpr)

# ==========================================================
# heatmap plot
# ==========================================================
def plot_heatmap(
    seq,
    ncpr,
    protein_name,
    window,
    fig_height,
    tick_interval,
    title_size,
    label_size,
    tick_size,
    font_weight
):

    fig = plt.figure(figsize=(12, fig_height))

    gs = gridspec.GridSpec(
        1, 2,
        width_ratios=[45, 2],
        wspace=0.05
    )

    ax = fig.add_subplot(gs[0])
    cax = fig.add_subplot(gs[1])

    im = ax.imshow(
        ncpr.reshape(1, -1),
        cmap="seismic",
        vmin=-1,
        vmax=1,
        aspect="auto"
    )

    # =========================
    # axis
    # =========================
    ax.set_yticks([])
    ax.set_xlabel("Residue Position", fontsize=label_size, fontweight=font_weight)

    ax.set_title(
        f"{protein_name} NCPR (window={window})",
        fontsize=title_size,
        fontweight=font_weight
    )

    ticks = np.arange(0, len(seq) + 1, tick_interval)
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks], fontsize=tick_size)

    # =========================
    # colorbar
    # =========================
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("Average Charge (NCPR)", fontsize=label_size, fontweight=font_weight)
    cbar.ax.tick_params(labelsize=tick_size)

    return fig

# ==========================================================
# MAIN UI
# ==========================================================
st.title("Protein NCPR Heatmap Tool")

protein_name = st.text_input("Protein name", "ProteinA")
seq_input = st.text_area("Amino acid sequence")

window = st.number_input("Sliding window size", 3, 101, 21, 2)

tick_interval = st.number_input("X-axis tick interval", 1, 200, 20)

fig_height = st.slider("Heatmap height", 1.0, 4.0, 1.6, 0.1)

# =========================
# font settings
# =========================
st.markdown("### Font settings")

title_size = st.slider("Title font size", 8, 20, 10)
label_size = st.slider("Axis label font size", 8, 20, 9)
tick_size = st.slider("Tick font size", 6, 18, 8)

font_weight = st.selectbox("Font weight", ["normal", "bold"])

# ==========================================================
# RUN
# ==========================================================
if st.button("Generate"):

    seq = seq_input.replace(" ", "").replace("\n", "").upper()

    if len(seq) == 0:
        st.error("Please input sequence")
    elif len(seq) < window:
        st.error("Sequence shorter than window")
    else:
        ncpr = compute_ncpr(seq, window)

        fig = plot_heatmap(
            seq,
            ncpr,
            protein_name,
            window,
            fig_height,
            tick_interval,
            title_size,
            label_size,
            tick_size,
            font_weight
        )

        st.pyplot(fig)

        fig.savefig("ncpr_heatmap.png", dpi=300, bbox_inches="tight")

        with open("ncpr_heatmap.png", "rb") as f:
            st.download_button(
                "Download PNG",
                f,
                file_name="ncpr_heatmap.png"
            )