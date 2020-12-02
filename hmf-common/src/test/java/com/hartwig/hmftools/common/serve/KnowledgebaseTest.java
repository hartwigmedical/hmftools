package com.hartwig.hmftools.common.serve;

import static org.junit.Assert.assertEquals;

import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import org.junit.Test;

public class KnowledgebaseTest {

    @Test
    public void canJoinSources() {
        Set<Knowledgebase> sources = Sets.newHashSet(Lists.newArrayList(Knowledgebase.VICC_CGI, Knowledgebase.VICC_CIVIC));
        assertEquals("vicc_cgi,vicc_civic", Knowledgebase.commaSeparatedSourceString(sources));

        Set<Knowledgebase> sources2 = Sets.newHashSet(Lists.newArrayList(Knowledgebase.VICC_CIVIC, Knowledgebase.VICC_CGI));
        assertEquals("vicc_cgi,vicc_civic", Knowledgebase.commaSeparatedSourceString(sources2));
    }

}