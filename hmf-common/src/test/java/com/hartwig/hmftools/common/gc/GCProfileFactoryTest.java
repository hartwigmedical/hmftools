package com.hartwig.hmftools.common.gc;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;

import com.google.common.collect.Multimap;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.chromosome.Chromosome;

import org.junit.Test;

public class GCProfileFactoryTest {

    private static final String BASE_PATH = Resources.getResource("gc").getPath() + File.separator;

    @Test
    public void canLoadNormalFile() throws IOException {
        final Multimap<Chromosome, GCProfile> gcContent = GCProfileFactory.loadGCContent(1000, BASE_PATH + "GC_profile.1000bp.cnp");
        assertEquals(100, gcContent.size());
    }
}
