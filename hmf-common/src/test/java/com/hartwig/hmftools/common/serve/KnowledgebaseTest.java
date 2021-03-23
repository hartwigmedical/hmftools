package com.hartwig.hmftools.common.serve;

import static org.junit.Assert.assertEquals;

import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import org.junit.Test;

public class KnowledgebaseTest {

    @Test
    public void canConvertSources() {
        Set<Knowledgebase> sources = Sets.newHashSet(Lists.newArrayList(Knowledgebase.VICC_CGI, Knowledgebase.VICC_CIVIC));
        String sourceString = Knowledgebase.toCommaSeparatedSourceString(sources);
        assertEquals("vicc_cgi,vicc_civic", sourceString);

        Set<Knowledgebase> convertedSources = Knowledgebase.fromCommaSeparatedTechnicalDisplayString(sourceString);
        assertEquals(sources, convertedSources);

        // Converting also sorts alphabetically:
        Set<Knowledgebase> sources2 = Sets.newHashSet(Lists.newArrayList(Knowledgebase.VICC_CIVIC, Knowledgebase.VICC_CGI));
        assertEquals("vicc_cgi,vicc_civic", Knowledgebase.toCommaSeparatedSourceString(sources2));
    }
}