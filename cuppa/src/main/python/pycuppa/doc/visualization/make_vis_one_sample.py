import importlib_resources as impresources
from cuppa.tests.mock_data import MockVisData

DOC_DIR = impresources.files("doc")

MockVisData.predictions.plot(
    sample_id="Breast_TNBC",
    plot_path=DOC_DIR/"visualization/cuppa_vis.png"
)

