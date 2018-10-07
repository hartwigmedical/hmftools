package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Optional;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CanonicalAnnotationTest {

    @Test
    public void favourDriverCatalogGenes() {
        final TranscriptAnnotation nonCanonicalDriverGene = new SimpleTranscriptAnnotation("ATP1A1", "ENST00000295598");
        final TranscriptAnnotation noDriverGene = new SimpleTranscriptAnnotation("AL136376.1", "ENST00000598661");
        final TranscriptAnnotation canonicalDriverGene = new SimpleTranscriptAnnotation("ATP1A1", "ENST00000537345");

        final CanonicalAnnotation victim = new CanonicalAnnotation();
        assertEquals(Optional.empty(), victim.pickCanonicalFavourDriverGene(Lists.newArrayList(nonCanonicalDriverGene)));

        Optional<TranscriptAnnotation> annotationSecond =
                victim.pickCanonicalFavourDriverGene(Lists.newArrayList(nonCanonicalDriverGene, noDriverGene));
        assertTrue(annotationSecond.isPresent());
        assertEquals(noDriverGene, annotationSecond.get());

        Optional<TranscriptAnnotation> annotationThird =
                victim.pickCanonicalFavourDriverGene(Lists.newArrayList(nonCanonicalDriverGene, noDriverGene, canonicalDriverGene));
        assertTrue(annotationThird.isPresent());
        assertEquals(canonicalDriverGene, annotationThird.get());
    }

    private static class SimpleTranscriptAnnotation implements TranscriptAnnotation {

        @NotNull
        private final String gene;
        @NotNull
        private final String transcript;

        SimpleTranscriptAnnotation(@NotNull final String gene, @NotNull final String transcript) {
            this.gene = gene;
            this.transcript = transcript;
        }

        @NotNull
        @Override
        public String gene() {
            return gene;
        }

        @NotNull
        @Override
        public String transcript() {
            return transcript;
        }
    }
}
