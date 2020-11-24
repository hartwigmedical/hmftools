package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.serve.classification.MutationType;

import org.junit.Test;

public class MutationTypeExtractorTest {

    @Test
    public void canDetermineMutationTypes() {
        MutationTypeExtractor mutationTypeExtractor = MutationTypeExtractor.buildProductionExtractor();
        assertEquals(MutationType.HOTSPOT, mutationTypeExtractor.extractType("BRAF", "V600E"));

        assertEquals(MutationType.UNKNOWN, mutationTypeExtractor.extractType("BRAF", "what is this?"));

        // If gene symbol is missing, we always classify as UNKNOWN.
        assertEquals(MutationType.UNKNOWN, mutationTypeExtractor.extractType(null, "V600E"));
    }
}