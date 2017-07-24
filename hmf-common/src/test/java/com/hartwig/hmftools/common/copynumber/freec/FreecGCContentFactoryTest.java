package com.hartwig.hmftools.common.copynumber.freec;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;

import com.google.common.collect.Multimap;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.purple.ratio.GCContent;

import org.junit.Test;

public class FreecGCContentFactoryTest {

    private static final String BASE_PATH = Resources.getResource("copynumber").getPath() + File.separator + "gccontent";

    @Test
    public void canLoadNormalFile() throws IOException, HartwigException {
        final Multimap<String, GCContent> gcContent = FreecGCContentFactory.loadGCContent(BASE_PATH);
        assertEquals(100, gcContent.size());
    }
}
