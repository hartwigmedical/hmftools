package com.hartwig.hmftools.serve.vicc;

import static org.junit.Assert.assertEquals;

import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import org.junit.Test;

public class ViccExtractorTestApplicationTest {

    @Test
    public void canJoinSources() {
        Set<String> sources = Sets.newTreeSet(Lists.newArrayList("SRC1", "SRC2", "SRC3"));
        assertEquals("SRC1,SRC2,SRC3", ViccExtractorTestApplication.buildSourcesString(sources));
    }
}