package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class ViccEventTypeExtractorTest {

    @Test
    public void canDetermineMutationTypes() {
        assertEquals(EventType.HOTSPOT, ViccEventTypeExtractor.extractType(createFeature("BRAF", "V600E")));

        assertEquals(EventType.UNKNOWN, ViccEventTypeExtractor.extractType(createFeature("BRAF", "what is this?")));

        // If gene symbol is missing, we always classify as UNKNOWN.
        assertEquals(EventType.UNKNOWN, ViccEventTypeExtractor.extractType(createFeature(null, "V600E")));
    }

    @NotNull
    private static Feature createFeature(@Nullable String gene, @NotNull String name) {
        return ImmutableFeature.builder().geneSymbol(gene).name(name).build();
    }
}