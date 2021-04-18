package com.hartwig.hmftools.common.genome.region;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.collect.SortedSetMultimap;
import com.google.common.io.Resources;

import org.junit.Test;

public class BEDFileLoaderTest {

    private static final String EXAMPLE_BED_FILE = Resources.getResource("bed/example.bed").getPath();

    @Test
    public void verifyStartIsOneBased() throws IOException {
        SortedSetMultimap<String, GenomeRegion> regions = BEDFileLoader.fromBedFile(EXAMPLE_BED_FILE);
        assertEquals(1, regions.get("1").first().start());
        assertEquals(1, regions.get("1").first().end());
    }
}
