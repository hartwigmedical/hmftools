package com.hartwig.hmftools.common.serve;

import static org.junit.Assert.assertEquals;

import java.util.Set;

import com.google.common.collect.Sets;

import org.junit.Test;

public class KnowledgebaseTest {

    @Test
    public void canConvertBackAndForth() {
        Set<Knowledgebase> sources = Sets.newHashSet(Knowledgebase.ACTIN, Knowledgebase.CKB);
        assertEquals(sources, Knowledgebase.fromCommaSeparatedSourceString(Knowledgebase.toCommaSeparatedSourceString(sources)));
    }

    @Test
    public void canConvertSourcesToString() {
        Set<Knowledgebase> sources = Sets.newHashSet(Knowledgebase.VICC_CGI, Knowledgebase.VICC_CIVIC);
        String sourceString = Knowledgebase.toCommaSeparatedSourceString(sources);
        assertEquals("VICC_CGI,VICC_CIVIC", sourceString);

        Set<Knowledgebase> convertedSources = Knowledgebase.fromCommaSeparatedSourceString(sourceString);
        assertEquals(sources, convertedSources);

        // Converting also sorts alphabetically:
        Set<Knowledgebase> sources2 = Sets.newHashSet(Knowledgebase.VICC_CIVIC, Knowledgebase.VICC_CGI);
        assertEquals("VICC_CGI,VICC_CIVIC", Knowledgebase.toCommaSeparatedSourceString(sources2));
    }
}