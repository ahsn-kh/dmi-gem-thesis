from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

EXPECTED_FIGURES = {
    "Figure_2_1.png",
    "Figure_2_2.png",
    "Figure_3_1.png",
    "Figure_3_2.png",
    "Figure_3_3.png",
    "Figure_3_4.png",
    "Figure_3_5.png",
    "Figure_3_6.png",
    "Figure_4_1.png",
    "Figure_A1.png",
    "Figure_A2.png",
    "Figure_A3.png",
}


def test_public_raw_inputs_present():
    assert (ROOT / "data/raw_public/GSE17576/GSE17576_series_matrix.txt.gz").exists()
    assert (ROOT / "data/raw_public/GSE305719/GSE305719_data_count.csv.gz").exists()


def test_processed_transcriptomes_present():
    for dataset in ("GSE17576", "GSE305719"):
        base = ROOT / "data/processed_transcriptomics" / dataset
        assert (base / "normalized_expression.csv.xz").exists()
        assert (base / "sample_annotation.csv.xz").exists()


def test_transcriptome_processing_code_present():
    assert (ROOT / "scripts/transcriptomics/prepare_transcriptome_datasets.py").exists()
    assert (ROOT / "scripts/transcriptomics/load_datasets.py").exists()


def test_current_thesis_figures_exact_set():
    actual = {p.name for p in (ROOT / "figures").glob("*.png")}
    assert actual == EXPECTED_FIGURES


def test_no_macos_metadata():
    assert not list(ROOT.rglob("__MACOSX"))
    assert not list(ROOT.rglob("._*"))


def test_dmi_signal_table_present():
    assert (ROOT / "data/dmi/table_voxel_met_conc.csv.xz").exists()


def test_exact_week9_dmi_target_table_present():
    assert (
        ROOT
        / "data/dmi"
        / "condition_summary_used__met1_modelFitC6_c6_aware_fva_sampling_atp_per_replicate_n50.csv.gz"
    ).exists()



def test_analysis_ready_dmi_input_present():
    assert (
        ROOT
        / "data"
        / "dmi"
        / "data_tissue_vals_model_fitting.joblib.xz"
    ).exists()


def test_dmi_extraction_and_fitting_code_present():
    assert (
        ROOT
        / "scripts"
        / "upstream_dmi"
        / "extraction"
        / "extract_met_conc.ipynb"
    ).exists()

    assert (
        ROOT
        / "scripts"
        / "upstream_dmi"
        / "fitting"
        / "model_fittingC6_improved.py"
    ).exists()
