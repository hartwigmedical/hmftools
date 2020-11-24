package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.serve.classification.MutationType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class MutationTypeExtractorTest {

    @Test
    public void canDetermineMutationTypes() {
        assertEquals(MutationType.HOTSPOT, MutationTypeExtractor.extractType(createFeature("BRAF", "V600E")));

        assertEquals(MutationType.UNKNOWN, MutationTypeExtractor.extractType(createFeature("BRAF", "what is this?")));

        // If gene symbol is missing, we always classify as UNKNOWN.
        assertEquals(MutationType.UNKNOWN, MutationTypeExtractor.extractType(createFeature(null, "V600E")));
    }

    @NotNull
    private static Feature createFeature(@Nullable String gene, @NotNull String name) {
        return ImmutableFeature.builder().geneSymbol(gene).name(name).build();
    }
}