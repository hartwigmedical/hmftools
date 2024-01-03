package com.hartwig.hmftools.orange.algo.cuppa;

import static com.hartwig.hmftools.common.cuppa.CategoryType.COMBINED;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.cuppa.CuppaDataFile;
import com.hartwig.hmftools.common.cuppa.DataTypes;
import com.hartwig.hmftools.common.cuppa.ResultType;
import com.hartwig.hmftools.common.cuppa2.CuppaPredictionEntry;
import com.hartwig.hmftools.common.cuppa2.CuppaPredictions;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.cuppa.CuppaPrediction;
import com.hartwig.hmftools.datamodel.cuppa.ImmutableCuppaPrediction;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CuppaDataFactoryTest
{
    private static final double EPSILON = 1.0E-10;

    private static final String CUPPA_VIS_DATA_TSV = Resources.getResource("test_run/cuppa/tumor_sample.cuppa.vis_data.tsv").getPath();

    @Test
    public void canExtractPredictionsFromCuppaV2Output() throws IOException
    {
        CuppaPredictions cuppaPredictions = CuppaPredictions.fromTsv(CUPPA_VIS_DATA_TSV);
        CuppaData cuppaData = CuppaDataFactory.create(cuppaPredictions, null);
        CuppaPrediction prediction = cuppaData.predictions().get(0);

        assertEquals(prediction.cuppaMajorVersion(), "v2");
        assertEquals(prediction.cancerType(), "Breast: Triple negative");
        assertEquals(prediction.likelihood(), 0.7648, EPSILON);

        assertEquals(prediction.genomicPositionClassifier(), 0.6103, EPSILON);
        assertEquals(prediction.snv96Classifier(), 0.8326, EPSILON);
        assertEquals(prediction.eventClassifier(), 0.2607, EPSILON);

        assertEquals(prediction.geneExpressionClassifier(), 0.5215, EPSILON);
        assertEquals(prediction.altSjClassifier(), 0.6721, EPSILON);
    }

    private static final String CUPPA_TEST_CSV = Resources.getResource("test_run/cuppa/tumor_sample.cup.data.csv").getPath();

    @Test
    public void canCreateCorrectDataForTestRun() throws IOException
    {
        CuppaData cuppaData = CuppaDataFactory.create(
                null,
                CuppaDataFile.read(CUPPA_TEST_CSV)
        );

        assertEquals(36, cuppaData.predictions().size());
        Map<String, CuppaPrediction> actualPredictionsByCancerType =
                cuppaData.predictions().stream().collect(Collectors.toMap(CuppaPrediction::cancerType, entry -> entry));

        for(CuppaPrediction expected : List.of(
                prediction("v0", "Melanoma", 0.996, 0.979, 0.990, 0.972),
                prediction("v0","Pancreas", 0.00013, 6.75e-28, 6.05e-05, 1.8e-05),
                prediction("v0","Lung: Non-small Cell", 0.00013, 4.83e-28, 2.91e-05, 0.00518)))
        {
            CuppaPrediction actual = actualPredictionsByCancerType.get(expected.cancerType());
            assertNotNull(actual);
            List<Function<CuppaPrediction, Double>> functionsToVerify = List.of(CuppaPrediction::likelihood,
                    CuppaPrediction::snv96Classifier,
                    CuppaPrediction::genomicPositionClassifier,
                    CuppaPrediction::eventClassifier,
                    CuppaPrediction::altSjClassifier,
                    CuppaPrediction::geneExpressionClassifier);
            for(Function<CuppaPrediction, Double> function : functionsToVerify)
            {
                assertCuppaPredictionField(expected, actual, function);
            }
        }
    }

    @Test
    public void respectOrderingOfCombinedDataTypes()
    {
        CuppaDataFile rna = create(DataTypes.DATA_TYPE_RNA_COMBINED, "rna");
        CuppaDataFile dna = create(DataTypes.DATA_TYPE_DNA_COMBINED, "dna");
        CuppaDataFile overall = create(DataTypes.DATA_TYPE_COMBINED, "overall");

        CuppaData fromRna = CuppaDataFactory.create(null, Lists.newArrayList(rna));
        assertEquals("rna", fromRna.predictions().get(0).cancerType());

        CuppaData fromRnaDna = CuppaDataFactory.create(null, Lists.newArrayList(rna, dna));
        assertEquals("dna", fromRnaDna.predictions().get(0).cancerType());

        CuppaData fromAll = CuppaDataFactory.create(null, Lists.newArrayList(rna, dna, overall));
        assertEquals("overall", fromAll.predictions().get(0).cancerType());
    }

    @Test
    public void doNotCrashOnMissingEntries()
    {
        CuppaPredictions cuppaPredictionsV2 = new CuppaPredictions(new ArrayList<>());
        assertNotNull(CuppaDataFactory.create(cuppaPredictionsV2,null));

        List<CuppaDataFile> cuppaPredictionsV1 = new ArrayList<>();
        assertNotNull(CuppaDataFactory.create(null, cuppaPredictionsV1));
    }

    @NotNull
    private static ImmutableCuppaPrediction prediction(
            @NotNull String cuppaMajorVersion,
            @NotNull String cancerType,
            double likelihood,
            double snvPairwiseClassifier,
            double genomicPositionClassifier,
            double eventClassifier
    ){
        return ImmutableCuppaPrediction.builder()
                .cuppaMajorVersion(cuppaMajorVersion)
                .cancerType(cancerType)
                .likelihood(likelihood)
                .snv96Classifier(snvPairwiseClassifier)
                .genomicPositionClassifier(genomicPositionClassifier)
                .eventClassifier(eventClassifier)
                .build();
    }

    private static void assertCuppaPredictionField(@NotNull CuppaPrediction expected, @NotNull CuppaPrediction actual,
            @NotNull Function<CuppaPrediction, Double> function)
    {
        Double expectedValue = function.apply(expected);
        Double actualValue = function.apply(actual);
        if(expectedValue == null)
        {
            assertNull(actualValue);
        }
        else
        {
            assertNotNull(actualValue);
            assertEquals(expectedValue, actualValue, EPSILON);
        }
    }

    @NotNull
    private static CuppaDataFile create(@NotNull String dataType, @NotNull String refCancerType)
    {
        Map<String, Double> cancerTypeValues = Maps.newHashMap();
        cancerTypeValues.put(refCancerType, 0D);

        return new CuppaDataFile(COMBINED, ResultType.CLASSIFIER, dataType, Strings.EMPTY, cancerTypeValues);
    }
}