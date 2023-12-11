import pandas as pd
from sklearn.compose import make_column_selector

from cuppa.constants import DEFAULT_FUSION_OVERRIDES_PATH
from cuppa.classifier.classifiers import SubClassifiers, MetaClassifiers, ClassifierLayers
from cuppa.misc.mock_data import MockTrainingData, MockProbsFromFitTransform


def get_expected_probs(clf_names: str | list[str]) -> pd.DataFrame:
    clf_names = pd.Series(clf_names)

    if len(clf_names)==1:
        return getattr(MockProbsFromFitTransform, clf_names[0])

    probs_dict = {}
    for clf_name in clf_names:
        probs_dict[clf_name] = getattr(MockProbsFromFitTransform, clf_name)

    probs_merged = []
    for clf_name, probs in probs_dict.items():
        probs = probs.copy()
        probs.columns = clf_name + "__" + probs.columns
        probs_merged.append(probs)

    probs_merged = pd.concat(probs_merged, axis=1)

    return probs_merged


class TestSubClassifiers:
    #from cuppa.classifier.test.test_classifiers import TestSubClassifiers
    #self = TestSubClassifiers

    X = MockTrainingData.X
    y = MockTrainingData.y

    def test_fit_transform_gen_pos(self):
        classifier = SubClassifiers.GenPosClassifier()
        probs = classifier.fit_transform(self.X[make_column_selector("^gen_pos")], self.y).round(3)
        probs_expected = get_expected_probs("gen_pos").round(3)
        assert probs.equals(probs_expected)

    def test_fit_transform_snv96(self):
        classifier = SubClassifiers.Snv96Classifier()
        probs = classifier.fit_transform(self.X[make_column_selector("^snv96")], self.y).round(3)
        probs_expected = MockProbsFromFitTransform.snv96.round(3)
        assert probs.equals(probs_expected)

    def test_fit_transform_event(self):
        classifier = SubClassifiers.EventClassifier()
        probs = classifier.fit_transform(self.X[make_column_selector("^event")], self.y).round(3)
        probs_expected = MockProbsFromFitTransform.event.round(3)
        assert probs.equals(probs_expected)

    def test_fit_transform_gene_exp(self):
        classifier = SubClassifiers.GeneExpClassifier()
        probs = classifier.fit_transform(self.X[make_column_selector("^gene_exp")], self.y).round(3)
        probs_expected = MockProbsFromFitTransform.gene_exp.round(3)
        assert probs.equals(probs_expected)

    def test_fit_transform_alt_sj(self):
        classifier = SubClassifiers.AltSjClassifier()
        probs = classifier.fit_transform(self.X[make_column_selector("^alt_sj")], self.y).round(3)
        probs_expected = MockProbsFromFitTransform.alt_sj.round(3)
        assert probs.equals(probs_expected)

    def test_fit_transform_sub_classifier_layer(self):
        classifier = ClassifierLayers.SubClassifierLayer()
        probs = classifier.fit_transform(self.X, self.y).round(3)
        probs_expected = get_expected_probs(["gen_pos", "snv96", "event", "gene_exp", "alt_sj"]).round(3)
        assert probs.equals(probs_expected)


class TestMetaClassifiers:
    # from cuppa.classifier.test.test_classifiers import TestMetaClassifiers
    # self = TestMetaClassifiers

    X = MockTrainingData.X
    y = MockTrainingData.y

    def test_fit_transform_dna_combined(self):
        classifier = MetaClassifiers.DnaCombinedClassifier(fusion_overrides_path=DEFAULT_FUSION_OVERRIDES_PATH)

        ## Set up overrides
        classifier["sex_filter"].sample_sexes = MockTrainingData.X[make_column_selector("is_male")].iloc[:,0].astype(bool)
        classifier["fusion_overrider"].sample_fusions = MockTrainingData.X[make_column_selector("fusion")]

        probs_sub_clf = get_expected_probs(["gen_pos", "snv96", "event"])

        probs = classifier.fit_transform(probs_sub_clf, MockTrainingData.y).round(3)
        probs_expected = get_expected_probs("dna_combined").round(3)

        assert probs.equals(probs_expected)

    def test_fit_transform_rna_combined(self):
        classifier = MetaClassifiers.RnaCombinedClassifier()

        ## Set up overrides
        classifier["sex_filter"].sample_sexes = MockTrainingData.X[make_column_selector("is_male")].iloc[:, 0].astype(bool)

        probs_sub_clf = get_expected_probs(["gene_exp", "alt_sj"])
        probs_sub_clf = probs_sub_clf.reindex(MockTrainingData.y.index)

        probs_expected = get_expected_probs("rna_combined").round(3)

        probs = classifier.fit_transform(probs_sub_clf, MockTrainingData.y).round(3)
        probs = probs.reindex(probs_expected.index).reindex(probs_expected.columns, axis=1)


        assert probs.equals(probs_expected)
