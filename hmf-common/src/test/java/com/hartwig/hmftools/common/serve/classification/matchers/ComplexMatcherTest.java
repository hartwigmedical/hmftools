package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ComplexMatcherTest {

    private static final Map<String, Set<String>> COMPLEX_EVENTS_PER_GENE = Maps.newHashMap();

    static {
        COMPLEX_EVENTS_PER_GENE.put("VHL", Sets.newHashSet("Splicing alteration (c.464-2A>G)"));
    }

    @Test
    public void canAssessWhetherEventIsComplexEvent() {
        EventMatcher matcher = new ComplexMatcher(testHotspotMatcher(), COMPLEX_EVENTS_PER_GENE);

        assertTrue(matcher.matches("VHL", "Splicing alteration (c.464-2A>G)"));
        assertTrue(matcher.matches("KRAS", "KRAS ."));
        assertTrue(matcher.matches("APC", "APC p.I1557*fs*1"));
        assertTrue(matcher.matches("BRCA1", "BRCA1 L631Qfs*4"));

        assertFalse(matcher.matches("APC", "APC p.I1557fs"));

        assertTrue(matcher.matches("ERBB2", "S310F/Y"));
        assertTrue(matcher.matches("ERBB2", "L698_S1037dup"));
        assertFalse(matcher.matches("EGFR", "Exon 19 deletion/insertion"));

        assertFalse(matcher.matches("BRCA1", "BRCA1 L631QFS"));
    }

    @NotNull
    private static HotspotMatcher testHotspotMatcher() {
        return new HotspotMatcher(event -> event, new FusionPairMatcher(Sets.newHashSet(), Sets.newHashSet(), Sets.newHashSet()));
    }
}