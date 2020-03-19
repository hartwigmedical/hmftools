package com.hartwig.hmftools.knowledgebasegenerator.sourceknowledgebase;

import static org.junit.Assert.*;

import org.junit.Test;

public class SourceTest {

    @Test
    public void canDetermineSource() {
        assertEquals(Source.ONCOKB, Source.sourceFromKnowledgebase("oncokb"));
        assertEquals(Source.CGI, Source.sourceFromKnowledgebase("cgi"));
        assertEquals(Source.CIVIC, Source.sourceFromKnowledgebase("civic"));
        assertEquals(Source.JAX, Source.sourceFromKnowledgebase("jax"));
        assertEquals(Source.JAX_TRIALS, Source.sourceFromKnowledgebase("jax_trials"));
        assertEquals(Source.BRCA, Source.sourceFromKnowledgebase("brca"));
        assertEquals(Source.SAGE, Source.sourceFromKnowledgebase("sage"));
        assertEquals(Source.PMKB, Source.sourceFromKnowledgebase("pmkb"));
        assertEquals(Source.MOLECULARMATCH, Source.sourceFromKnowledgebase("molecularmatch"));
        assertEquals(Source.MOLECULARMATCH_TRIALS, Source.sourceFromKnowledgebase("molecularmatch_trials"));
    }
}