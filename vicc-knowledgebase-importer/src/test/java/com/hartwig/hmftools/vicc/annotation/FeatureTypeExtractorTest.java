package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FeatureTypeExtractorTest {

    @Test
    public void canDetermineFeatureTypes() {
        assertEquals(FeatureType.HOTSPOT, FeatureTypeExtractor.extractType(createFeatureWithNameAndGene("V600E", "BRAF")));

        assertEquals(FeatureType.UNKNOWN, FeatureTypeExtractor.extractType(createFeatureWithNameAndGene("what is this?", "BRAF")));
    }

    @NotNull
    private static Feature createFeatureWithNameAndGene(@NotNull String name, @NotNull String gene) {
        return ImmutableFeature.builder().name(name).geneSymbol(gene).build();
    }

}