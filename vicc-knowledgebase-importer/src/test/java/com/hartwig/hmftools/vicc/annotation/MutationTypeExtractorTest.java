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
        assertEquals(MutationType.HOTSPOT, MutationTypeExtractor.extractType(createFeatureWithGeneAndName("BRAF", "V600E")));

        assertEquals(MutationType.UNKNOWN, MutationTypeExtractor.extractType(createFeatureWithGeneAndName("BRAF", "what is this?")));

        assertEquals(MutationType.UNKNOWN, MutationTypeExtractor.extractType(createFeatureWithGeneAndName(null, "V600E")));
    }

    @NotNull
    private static Feature createFeatureWithGeneAndName(@Nullable String gene, @NotNull String name) {
        return ImmutableFeature.builder().name(name).geneSymbol(gene).build();
    }
}