package com.hartwig.hmftools.orange.algo.cuppa;

import static com.hartwig.hmftools.common.cuppa.CategoryType.COMBINED;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.cuppa.CuppaDataFile;
import com.hartwig.hmftools.common.cuppa.DataTypes;
import com.hartwig.hmftools.common.cuppa.ResultType;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CuppaPredictionFactoryTest {

    private static final String CUPPA_TEST_CSV = Resources.getResource("test_run/cuppa/tumor_sample.cup.data.csv").getPath();

    private static final double EPSILON = 1.0E-10;

    @Test
    public void canCreateCorrectDataForTestRun() throws IOException {
        CuppaData cuppaData = CuppaDataFactory.create(CuppaDataFile.read(CUPPA_TEST_CSV));

        assertEquals(36, cuppaData.predictions().size());
        assertEquals(0.996, findByCancerType(cuppaData.predictions(), "Melanoma").likelihood(), EPSILON);
        assertEquals(0.000131, findByCancerType(cuppaData.predictions(), "Pancreas").likelihood(), EPSILON);
        assertEquals(0.000146, findByCancerType(cuppaData.predictions(), "Lung: Non-small Cell").likelihood(), EPSILON);

        assertEquals(3, cuppaData.simpleDups32To200B());
        assertEquals(8, cuppaData.maxComplexSize());
        assertEquals(0, cuppaData.telomericSGLs());
        assertEquals(4, cuppaData.LINECount());
    }

    @NotNull
    private static CuppaPrediction findByCancerType(@NotNull List<CuppaPrediction> predictions, @NotNull String cancerType) {
        for (CuppaPrediction prediction : predictions) {
            if (prediction.cancerType().equals(cancerType)) {
                return prediction;
            }
        }

        throw new IllegalStateException("Could not find prediction for cancer type '" + cancerType + "'");
    }

    @Test
    public void respectOrderingOfCombinedDataTypes() {
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
    public void doNotCrashOnMissingEntries() {
        assertNotNull(CuppaDataFactory.create(Lists.newArrayList()));
    }

    @NotNull
    private static CuppaDataFile create(@NotNull String dataType, @NotNull String refCancerType) {
        Map<String, Double> cancerTypeValues = Maps.newHashMap();
        cancerTypeValues.put(refCancerType, 0D);

        return new CuppaDataFile(COMBINED, ResultType.CLASSIFIER, dataType, Strings.EMPTY, cancerTypeValues);
    }
}