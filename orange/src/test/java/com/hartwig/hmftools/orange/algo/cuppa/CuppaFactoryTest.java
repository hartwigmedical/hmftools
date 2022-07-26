package com.hartwig.hmftools.orange.algo.cuppa;

import static com.hartwig.hmftools.common.cuppa.CategoryType.COMBINED;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cuppa.CuppaEntry;
import com.hartwig.hmftools.common.cuppa.DataTypes;
import com.hartwig.hmftools.common.cuppa.ImmutableCuppaEntry;
import com.hartwig.hmftools.common.cuppa.ResultType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CuppaFactoryTest {

    @Test
    public void respectOrderingOfCombinedDataTypes() {
        CuppaEntry rna = combined().dataType(DataTypes.DATA_TYPE_RNA_COMBINED).refCancerType("rna").build();
        CuppaEntry dna = combined().dataType(DataTypes.DATA_TYPE_DNA_COMBINED).refCancerType("dna").build();
        CuppaEntry overall = combined().dataType(DataTypes.DATA_TYPE_COMBINED).refCancerType("overall").build();

        CuppaData fromRna = CuppaFactory.create(Lists.newArrayList(rna));
        assertEquals("rna", fromRna.predictions().get(0).cancerType());

        CuppaData fromRnaDna = CuppaFactory.create(Lists.newArrayList(rna, dna));
        assertEquals("dna", fromRnaDna.predictions().get(0).cancerType());

        CuppaData fromAll = CuppaFactory.create(Lists.newArrayList(rna, dna, overall));
        assertEquals("overall", fromAll.predictions().get(0).cancerType());
    }

    @Test
    public void doNotCrashOnMissingEntries() {
        assertNotNull(CuppaFactory.create(Lists.newArrayList()));
    }

    @NotNull
    private static ImmutableCuppaEntry.Builder combined() {
        return ImmutableCuppaEntry.builder().category(COMBINED).resultType(ResultType.CLASSIFIER).value("").refValue(0);
    }
}