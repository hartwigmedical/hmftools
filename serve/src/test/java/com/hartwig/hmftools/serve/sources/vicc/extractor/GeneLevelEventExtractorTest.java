package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static org.junit.Assert.*;

import com.hartwig.hmftools.serve.actionability.gene.GeneLevelEvent;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;

import org.junit.Test;

public class GeneLevelEventExtractorTest {

    @Test
    public void canExtractGeneLevelEvent() {
        Feature feature = ImmutableFeature.builder().geneSymbol("a").name("a").build();

        assertEquals(GeneLevelEventExtractor.extractGeneLevelEvent("promiscuousFusion", feature), GeneLevelEvent.FUSION);
        assertEquals(GeneLevelEventExtractor.extractGeneLevelEvent("geneLevel", feature), GeneLevelEvent.ACTIVATION);
        assertEquals(GeneLevelEventExtractor.extractGeneLevelEvent("ccc", feature), GeneLevelEvent.UNKONWN);
        assertNotEquals(GeneLevelEventExtractor.extractGeneLevelEvent("aa", feature), GeneLevelEvent.INACTIVATION);

    }

}