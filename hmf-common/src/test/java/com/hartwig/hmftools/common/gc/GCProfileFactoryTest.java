package com.hartwig.hmftools.common.gc;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;

import com.google.common.collect.Multimap;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HartwigException;

import org.junit.Test;

public class GCProfileFactoryTest {

    private static final String BASE_PATH = Resources.getResource("gcprofile").getPath() + File.separator;

    @Test
    public void canLoadNormalFile() throws IOException, HartwigException {
        final Multimap<String, GCProfile> gcContent = GCProfileFactory.loadGCContent(BASE_PATH + "GC_profile.1000bp.cnp");
        assertEquals(100, gcContent.size());
    }
}
