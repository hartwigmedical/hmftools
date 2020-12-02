package com.hartwig.hmftools.serve.hotspot;

import static org.junit.Assert.assertEquals;

import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.ServeTestFactory;

import org.junit.Test;

public class HotspotVCFTest {

    @Test
    public void canDetermineUniqueSources() {
        KnownHotspot hotspot1 = ServeTestFactory.createTestKnownHotspotForSource(Knowledgebase.CIVIC);
        KnownHotspot hotspot2 = ServeTestFactory.createTestKnownHotspotForSource(Knowledgebase.CGI);
        KnownHotspot hotspot3 = ServeTestFactory.createTestKnownHotspotForSource(Knowledgebase.ONCOKB);
        KnownHotspot hotspot4 = ServeTestFactory.createTestKnownHotspotForSource(Knowledgebase.CIVIC);

        assertEquals("cgi,civic,oncokb", HotspotVCF.uniqueSourcesString(Lists.newArrayList(hotspot1, hotspot2, hotspot3, hotspot4)));
    }

    @Test
    public void canJoinSources() {
        Set<Knowledgebase> sources = Sets.newTreeSet(Lists.newArrayList(Knowledgebase.CGI, Knowledgebase.CIVIC));
        assertEquals("cgi,civic", HotspotVCF.toSourceString(sources));

        Set<Knowledgebase> sources2 = Sets.newTreeSet(Lists.newArrayList(Knowledgebase.CIVIC, Knowledgebase.CGI));
        assertEquals("cgi,civic", HotspotVCF.toSourceString(sources2));
    }

}