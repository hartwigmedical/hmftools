import importlib
from cuppa.misc.mock_data import MockCvOutput

DOC_DIR = importlib.resources.files("doc")

MockCvOutput.predictions_for_vis.plot(
    sample_id="Breast_TNBC",
    vis_data_out_path=DOC_DIR/"visualization/cuppa_vis_data.tsv",
    plot_path=DOC_DIR/"visualization/cuppa_vis.png"
)

