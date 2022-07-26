package com.hartwig.hmftools.orange.algo.cuppa;

import static com.hartwig.hmftools.common.cuppa.CategoryType.COMBINED;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.cuppa.CuppaDataFile;
import com.hartwig.hmftools.common.cuppa.DataTypes;
import com.hartwig.hmftools.common.cuppa.ResultType;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CuppaDataFactoryTest {

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