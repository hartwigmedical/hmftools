from cuppa.sample_data.cuppa_features import CuppaFeaturesLoader
from cuppa.classifier.cuppa_classifier import CuppaClassifier
from cuppa.constants import DEFAULT_CUPPA_CLASSIFIER_PATH
from cuppa.visualization.visualization import CuppaVisDataBuilder, CuppaVisPlotter

import importlib_resources as impresources
INPUT_PATH = impresources.files("cuppa") / "resources/mock_data/input_data/new_format/COLO829v003T.cuppa_data.tsv.gz"
PLOT_PATH = impresources.files("doc") / "visualization/COLO829v003T.cuppa.vis.png"


loader = CuppaFeaturesLoader(str(INPUT_PATH), sample_id="COLO829v003T")
features = loader.load()

classifier = CuppaClassifier.from_file(DEFAULT_CUPPA_CLASSIFIER_PATH)
features = classifier.fill_missing_cols(features)

predictions = classifier.predict(features)

builder = CuppaVisDataBuilder(predictions)
vis_data = builder.build()

plotter = CuppaVisPlotter(vis_data, plot_path=str(PLOT_PATH))
plotter.plot()
