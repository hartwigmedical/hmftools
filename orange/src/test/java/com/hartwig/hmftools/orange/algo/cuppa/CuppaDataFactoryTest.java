package com.hartwig.hmftools.orange.algo.cuppa;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.cuppa2.CuppaPredictions;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.cuppa.CuppaPrediction;
import com.hartwig.hmftools.datamodel.cuppa.ImmutableCuppaPrediction;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class CuppaDataFactoryTest
{
    private static final double EPSILON = 1.0E-10;

    private static final String CUPPA_VIS_DATA_WITH_RNA_TSV = Resources.getResource("cuppa/with_rna.vis_data.tsv").getPath();
    private static final String CUPPA_VIS_DATA_WITHOUT_RNA_TSV =
            Resources.getResource("test_run/cuppa/tumor_sample.cuppa.vis_data.tsv").getPath();

    @Test
    public void canExtractPredictionsFromCuppaV2IncludingRna() throws Exception
    {
        CuppaPrediction expectedPredictionMelanoma = ImmutableCuppaPrediction.builder()
                .cancerType("Skin: Melanoma")
                .likelihood(0.9843)
                .genomicPositionClassifier(0.9643)
                .snvPairwiseClassifier(0.9632)
                .featureClassifier(0.9962)
                .expressionPairwiseClassifier(0.983)
                .altSjCohortClassifier(0.9691)
                .build();
        CuppaPrediction expectedPredictionSkinOther = ImmutableCuppaPrediction.builder()
                .cancerType("Skin: Other")
                .likelihood(0.0413136)
                .genomicPositionClassifier(0.0309102)
                .snvPairwiseClassifier(0.0040984)
                .featureClassifier(0.0377109)
                .expressionPairwiseClassifier(0.0106618)
                .altSjCohortClassifier(0.0291157)
                .build();
        CuppaPrediction expectedPredictionProstate = ImmutableCuppaPrediction.builder()
                .cancerType("Prostate")
                .likelihood(0.0325366)
                .genomicPositionClassifier(0.0229236)
                .snvPairwiseClassifier(0.0008691)
                .featureClassifier(0.0169241)
                .expressionPairwiseClassifier(0.041203)
                .altSjCohortClassifier(3.4526e-07)
                .build();

        Map<String, CuppaPrediction> expectedPredictionsByCancerType = new HashMap<>();
        expectedPredictionsByCancerType.put("Skin: Melanoma", expectedPredictionMelanoma);
        expectedPredictionsByCancerType.put("Skin: Other", expectedPredictionSkinOther);
        expectedPredictionsByCancerType.put("Prostate", expectedPredictionProstate);

        assertCuppaPredictions(CUPPA_VIS_DATA_WITH_RNA_TSV, expectedPredictionsByCancerType);
    }

    @Test
    public void canExtractPredictionsFromCuppaV2ExcludingRna() throws Exception
    {
        CuppaPrediction expectedPredictionMelanoma = ImmutableCuppaPrediction.builder()
                .cancerType("Skin: Melanoma")
                .likelihood(0.9999671059707099)
                .genomicPositionClassifier(0.9997492462021321)
                .snvPairwiseClassifier(0.9999879275632831)
                .featureClassifier(0.8507700198967394)
                .build();
        CuppaPrediction expectedPredictionSkinOther = ImmutableCuppaPrediction.builder()
                .cancerType("Skin: Other")
                .likelihood(0.0)
                .genomicPositionClassifier(3.579993902723518e-05)
                .snvPairwiseClassifier(8.659757942846464e-06)
                .featureClassifier(0.001662394830480488)
                .build();
        CuppaPrediction expectedPredictionProstate = ImmutableCuppaPrediction.builder()
                .cancerType("Prostate")
                .likelihood(0.0)
                .genomicPositionClassifier(3.1091986550511897e-06)
                .snvPairwiseClassifier(9.922426122823756e-16)
                .featureClassifier(0.0001238559647132041)
                .build();

        Map<String, CuppaPrediction> expectedPredictionsByCancerType = new HashMap<>();
        expectedPredictionsByCancerType.put("Skin: Melanoma", expectedPredictionMelanoma);
        expectedPredictionsByCancerType.put("Skin: Other", expectedPredictionSkinOther);
        expectedPredictionsByCancerType.put("Prostate", expectedPredictionProstate);

        assertCuppaPredictions(CUPPA_VIS_DATA_WITHOUT_RNA_TSV, expectedPredictionsByCancerType);
    }

    @Test
    public void canGetCorrectSvFeatureValueFromFileWithRna() throws Exception
    {
        CuppaPredictions cuppaPredictions = CuppaPredictions.fromTsv(CUPPA_VIS_DATA_WITH_RNA_TSV);
        int featureValue = CuppaDataFactory.getSvFeatureValue(cuppaPredictions, "sv.MAX_COMPLEX_SIZE");
        int expectedFeatureValue = 1000;
        assertEquals(expectedFeatureValue, featureValue);
    }

    @Test
    public void canGetCorrectSvFeatureValueFromFileWithoutRna() throws Exception
    {
        CuppaPredictions cuppaPredictions = CuppaPredictions.fromTsv(CUPPA_VIS_DATA_WITHOUT_RNA_TSV);
        int featureValue = CuppaDataFactory.getSvFeatureValue(cuppaPredictions, "sv.MAX_COMPLEX_SIZE");
        int expectedFeatureValue = 8;
        assertEquals(expectedFeatureValue, featureValue);
    }

    private static void assertCuppaPredictions(@NotNull String inputFileName,
            @NotNull Map<String, CuppaPrediction> expectedPredictionsByCancerType) throws Exception
    {
        CuppaPredictions cuppaPredictions = CuppaPredictions.fromTsv(inputFileName);
        CuppaData cuppaData = CuppaDataFactory.create(cuppaPredictions);
        List<CuppaPrediction> predictionEntries = cuppaData.predictions();

        assertEquals(40, predictionEntries.size());
        assertEquals("Skin: Melanoma", cuppaData.bestPrediction().cancerType());

        Map<String, CuppaPrediction> actualPredictionsByCancerType = new HashMap<>();
        for(CuppaPrediction prediction : predictionEntries)
        {
            actualPredictionsByCancerType.put(prediction.cancerType(), prediction);
        }

        for(String cancerType : expectedPredictionsByCancerType.keySet())
        {
            CuppaPrediction expectedPrediction = expectedPredictionsByCancerType.get(cancerType);
            CuppaPrediction actualPrediction = actualPredictionsByCancerType.get(cancerType);

            assertEquals(expectedPrediction.likelihood(), actualPrediction.likelihood(), EPSILON);
            assertEqualsNullableDouble(expectedPrediction.genomicPositionClassifier(), actualPrediction.genomicPositionClassifier());
            assertEqualsNullableDouble(expectedPrediction.snvPairwiseClassifier(), actualPrediction.snvPairwiseClassifier());
            assertEqualsNullableDouble(expectedPrediction.featureClassifier(), actualPrediction.featureClassifier());
            assertEqualsNullableDouble(expectedPrediction.expressionPairwiseClassifier(), actualPrediction.expressionPairwiseClassifier());
            assertEqualsNullableDouble(expectedPrediction.altSjCohortClassifier(), actualPrediction.altSjCohortClassifier());
        }
    }

    private static void assertEqualsNullableDouble(@Nullable Double expected, @Nullable Double actual)
    {
        if(expected == null)
        {
            assertNull(actual);
        }
        else
        {
            assertNotNull(actual);
            assertEquals(expected, actual, EPSILON);
        }
    }
}