package com.hartwig.hmftools.orange.algo.cuppa;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertThrows;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.cuppa.ClassifierName;
import com.hartwig.hmftools.common.cuppa.CuppaPredictions;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.cuppa.CuppaMode;
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
    public void canExtractPredictionsFromCuppaIncludingRna() throws Exception
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
        expectedPredictionsByCancerType.put(expectedPredictionMelanoma.cancerType(), expectedPredictionMelanoma);
        expectedPredictionsByCancerType.put(expectedPredictionSkinOther.cancerType(), expectedPredictionSkinOther);
        expectedPredictionsByCancerType.put(expectedPredictionProstate.cancerType(), expectedPredictionProstate);
        expectedPredictionsByCancerType.put("Bone/Soft tissue: Cartilaginous neoplasm", null);

        assertCuppaPredictions(CUPPA_VIS_DATA_WITH_RNA_TSV, 33, expectedPredictionsByCancerType);
    }

    @Test
    public void canExtractPredictionsFromCuppaExcludingRna() throws Exception
    {
        CuppaPrediction expectedPredictionMelanoma = ImmutableCuppaPrediction.builder()
                .cancerType("Skin: Melanoma")
                .likelihood(0.9999710902284965)
                .genomicPositionClassifier(0.9997621099072876)
                .snvPairwiseClassifier(0.9999910153956342)
                .featureClassifier(0.8991388522282147)
                .build();
        CuppaPrediction expectedPredictionSkinOther = ImmutableCuppaPrediction.builder()
                .cancerType("Skin: Other")
                .likelihood(0.0)
                .genomicPositionClassifier(3.159853588544685E-5)
                .snvPairwiseClassifier(5.697607828154083E-6)
                .featureClassifier(5.889241560947043E-4)
                .build();
        CuppaPrediction expectedPredictionProstate = ImmutableCuppaPrediction.builder()
                .cancerType("Prostate")
                .likelihood(0.0)
                .genomicPositionClassifier(2.8698904215136434E-6)
                .snvPairwiseClassifier(9.922426122823756e-16)
                .featureClassifier(3.744164257376693E-5)
                .build();
        CuppaPrediction expectedPredictionCartilaginousNeoplasm = ImmutableCuppaPrediction.builder()
                .cancerType("Bone/Soft tissue: Cartilaginous neoplasm")
                .likelihood(2.8909771503514605E-5)
                .genomicPositionClassifier(4.107161442240924e-10)
                .snvPairwiseClassifier(4.414947168645547e-16)
                .featureClassifier(7.756820351850527E-10)
                .build();

        Map<String, CuppaPrediction> expectedPredictionsByCancerType = new HashMap<>();
        expectedPredictionsByCancerType.put(expectedPredictionMelanoma.cancerType(), expectedPredictionMelanoma);
        expectedPredictionsByCancerType.put(expectedPredictionSkinOther.cancerType(), expectedPredictionSkinOther);
        expectedPredictionsByCancerType.put(expectedPredictionProstate.cancerType(), expectedPredictionProstate);
        expectedPredictionsByCancerType.put(expectedPredictionCartilaginousNeoplasm.cancerType(), expectedPredictionCartilaginousNeoplasm);

        assertCuppaPredictions(CUPPA_VIS_DATA_WITHOUT_RNA_TSV, 40, expectedPredictionsByCancerType);
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

    @Test
    public void canAssignCuppaModeFromFileWithoutRna() throws Exception
    {
        CuppaPredictions cuppaPredictions = CuppaPredictions.fromTsv(CUPPA_VIS_DATA_WITHOUT_RNA_TSV);
        CuppaMode mode = CuppaDataFactory.getCuppaMode(cuppaPredictions.MainCombinedClassifierName);
        CuppaMode expectedMode = CuppaMode.WGS;
        assertEquals(expectedMode, mode);
    }

    @Test
    public void canAssignCuppaModeFromFileWithRna() throws Exception
    {
        CuppaPredictions cuppaPredictions = CuppaPredictions.fromTsv(CUPPA_VIS_DATA_WITH_RNA_TSV);
        CuppaMode mode = CuppaDataFactory.getCuppaMode(cuppaPredictions.MainCombinedClassifierName);
        CuppaMode expectedMode = CuppaMode.WGTS;
        assertEquals(expectedMode, mode);
    }

    @Test
    public void canThrowExceptionForInvalidMainCombinedClassifierName()
    {
        assertThrows(IllegalArgumentException.class, ()-> CuppaDataFactory.getCuppaMode(ClassifierName.ALT_SJ));
    }

    private static void assertCuppaPredictions(@NotNull String inputFileName,
            int expectedPredictionSize,
            @NotNull Map<String, CuppaPrediction> expectedPredictionsByCancerType) throws Exception
    {
        CuppaData cuppaData = CuppaDataFactory.create(inputFileName);
        List<CuppaPrediction> predictionEntries = cuppaData.predictions();

        assertEquals(expectedPredictionSize, predictionEntries.size());
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
            if(expectedPrediction == null)
            {
                assertNull(actualPrediction);
            }
            else
            {
                assertEquals(expectedPrediction.likelihood(), actualPrediction.likelihood(), EPSILON);
                assertEqualsNullableDouble(expectedPrediction.genomicPositionClassifier(), actualPrediction.genomicPositionClassifier());
                assertEqualsNullableDouble(expectedPrediction.snvPairwiseClassifier(), actualPrediction.snvPairwiseClassifier());
                assertEqualsNullableDouble(expectedPrediction.featureClassifier(), actualPrediction.featureClassifier());
                assertEqualsNullableDouble(expectedPrediction.expressionPairwiseClassifier(), actualPrediction.expressionPairwiseClassifier());
                assertEqualsNullableDouble(expectedPrediction.altSjCohortClassifier(), actualPrediction.altSjCohortClassifier());
            }
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