package com.hartwig.hmftools.serve.hotspot.tools;

import static org.junit.Assert.assertEquals;

import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.junit.Test;

public class ServeHotspotGeneratorTest {

    @Test
    public void canJoinSources() {
        Set<Knowledgebase> sources = Sets.newTreeSet(Lists.newArrayList(Knowledgebase.CGI, Knowledgebase.CIVIC));
        assertEquals("cgi,civic", ServeHotspotGenerator.buildSourcesString(sources));
    }
}