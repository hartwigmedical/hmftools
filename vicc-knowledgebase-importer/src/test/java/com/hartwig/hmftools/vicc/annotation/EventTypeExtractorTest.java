package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class EventTypeExtractorTest {

    @Test
    public void canDetermineEventTypes() {
        assertEquals(EventType.HOTSPOT, EventTypeExtractor.extractType(createFeatureWithGeneAndName("BRAF", "V600E")));

        assertEquals(EventType.UNKNOWN, EventTypeExtractor.extractType(createFeatureWithGeneAndName("BRAF", "what is this?")));

        assertEquals(EventType.UNKNOWN, EventTypeExtractor.extractType(createFeatureWithGeneAndName(null, "V600E")));
    }

    @NotNull
    private static Feature createFeatureWithGeneAndName(@Nullable String gene, @NotNull String name) {
        return ImmutableFeature.builder().name(name).geneSymbol(gene).build();
    }
}