package com.hartwig.hmftools.common.genome.gc;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.collect.Multimap;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

import org.junit.Test;

public class GCProfileFactoryTest {

    private static final String GC_PROFILE_PATH = Resources.getResource("gc/GC_profile.1000bp.cnp").getPath();

    @Test
    public void canLoadNormalFile() throws IOException {
        final Multimap<Chromosome, GCProfile> gcContent = GCProfileFactory.loadGCContent(1000, GC_PROFILE_PATH);
        assertEquals(100, gcContent.size());
    }
}
