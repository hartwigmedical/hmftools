package com.hartwig.hmftools.orange.algo.cuppa;

import static com.hartwig.hmftools.common.cuppa.CategoryType.COMBINED;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.IOException;
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
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.cuppa.CuppaPrediction;
import com.hartwig.hmftools.datamodel.cuppa.ImmutableCuppaPrediction;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CuppaDataFactoryTest
{
    private static final String CUPPA_TEST_CSV = Resources.getResource("test_run/cuppa/tumor_sample.cup.data.csv").getPath();

    private static final double EPSILON = 1.0E-10;

    @Test
    public void canCreateCorrectDataForTestRun() throws IOException
    {
        CuppaData cuppaData = CuppaDataFactory.create(CuppaDataFile.read(CUPPA_TEST_CSV));

        assertEquals(36, cuppaData.predictions().size());
        Map<String, CuppaPrediction> actualPredictionsByCancerType =
                cuppaData.predictions().stream().collect(Collectors.toMap(CuppaPrediction::cancerType, entry -> entry));

        for(CuppaPrediction expected : List.of(prediction("Melanoma", 0.996, 0.979, 0.990, 0.972),
                prediction("Pancreas", 0.00013, 6.75e-28, 6.05e-05, 1.8e-05),
                prediction("Lung: Non-small Cell", 0.00013, 4.83e-28, 2.91e-05, 0.00518)))
        {
            CuppaPrediction actual = actualPredictionsByCancerType.get(expected.cancerType());
            assertNotNull(actual);
            List<Function<CuppaPrediction, Double>> functionsToVerify = List.of(CuppaPrediction::likelihood,
                    CuppaPrediction::snvPairwiseClassifier,
                    CuppaPrediction::genomicPositionClassifier,
                    CuppaPrediction::featureClassifier,
                    CuppaPrediction::altSjCohortClassifier,
                    CuppaPrediction::expressionPairwiseClassifier);
            for(Function<CuppaPrediction, Double> function : functionsToVerify)
            {
                assertCuppaPredictionField(expected, actual, function);
            }
        }

        assertEquals(3, cuppaData.simpleDups32To200B());
        assertEquals(8, cuppaData.maxComplexSize());
        assertEquals(0, cuppaData.telomericSGLs());
        assertEquals(3, cuppaData.lineCount());
    }

    @Test
    public void respectOrderingOfCombinedDataTypes()
    {
        CuppaDataFile rna = create(DataTypes.DATA_TYPE_RNA_COMBINED, "rna");
        CuppaDataFile dna = create(DataTypes.DATA_TYPE_DNA_COMBINED, "dna");
        CuppaDataFile overall = create(DataTypes.DATA_TYPE_COMBINED, "overall");

        CuppaData fromRna = CuppaDataFactory.create(Lists.newArrayList(rna));
        assertEquals("rna", fromRna.predictions().get(0).cancerType());

        CuppaData fromRnaDna = CuppaDataFactory.create(Lists.newArrayList(rna, dna));
        assertEquals("dna", fromRnaDna.predictions().get(0).cancerType());

        CuppaData fromAll = CuppaDataFactory.create(Lists.newArrayList(rna, dna, overall));
        assertEquals("overall", fromAll.predictions().get(0).cancerType());
    }

    @Test
    public void doNotCrashOnMissingEntries()
    {
        assertNotNull(CuppaDataFactory.create(Lists.newArrayList()));
    }

    @NotNull
    private static ImmutableCuppaPrediction prediction(@NotNull String cancerType, double likelihood, double snvPairwiseClassifier,
            double genomicPositionClassifier, double featureClassifier)
    {
        return ImmutableCuppaPrediction.builder()
                .cancerType(cancerType)
                .likelihood(likelihood)
                .snvPairwiseClassifier(snvPairwiseClassifier)
                .genomicPositionClassifier(genomicPositionClassifier)
                .featureClassifier(featureClassifier)
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