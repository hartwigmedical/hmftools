import importlib_resources as impresources
from cuppa.tests.mock_data import MockCvOutput

DOC_DIR = impresources.files("doc")

MockCvOutput.predictions_for_vis.plot(
    sample_id="Breast_TNBC",
    plot_path=DOC_DIR/"visualization/cuppa_vis.png"
)

