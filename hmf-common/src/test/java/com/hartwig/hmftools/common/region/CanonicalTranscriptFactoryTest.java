package com.hartwig.hmftools.common.region;

import static org.junit.Assert.assertEquals;

import java.io.InputStream;

import com.hartwig.hmftools.common.genepanel.HmfGenomeFileLoader;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CanonicalTranscriptFactoryTest {

    @Test
    public void testKRAS() {
        final HmfTranscriptRegion kras = select("KRAS.tsv");
        final CanonicalTranscript transcript = CanonicalTranscriptFactory.create(kras);
        assertTranscript(transcript, 6, 4, 1119, 189);
        assertEquals(25368375, transcript.codingStart());
        assertEquals(25398318, transcript.codingEnd());
    }

    @Test
    public void testPTEN() {
        final HmfTranscriptRegion region = select("PTEN.tsv");
        final CanonicalTranscript transcript = CanonicalTranscriptFactory.create(region);
        assertTranscript(transcript, 9, 9, 9027, 403);
    }

    @Test
    public void testWASH7P() {
        final HmfTranscriptRegion region = select("WASH7P.tsv");
        final CanonicalTranscript transcript = CanonicalTranscriptFactory.create(region);
        assertTranscript(transcript, 12, 0, 1783, 0);
    }

    private static void assertTranscript(@NotNull final CanonicalTranscript victim, int exons, int codingExons, long transcriptLength,
            int codons) {
        assertEquals(exons, victim.exons());
        assertEquals(codingExons, victim.codingExons());
        assertEquals(transcriptLength, victim.exonBases());
        assertEquals(3 * codons, victim.codingBases());
    }

    @NotNull
    private static HmfTranscriptRegion select(@NotNull final String file) {
        final InputStream inputStream = CanonicalTranscriptFactoryTest.class.getResourceAsStream("/gene/" + file);
        return HmfGenomeFileLoader.fromInputStream(inputStream).iterator().next();
    }
}
