package com.hartwig.hmftools.common.region.hmfslicer;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.SortedSet;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;

import org.junit.Test;

public class HmfGenomeFileLoaderTest {
    private static final String BASE_PATH = Resources.getResource("gene").getPath() + File.separator;
    private static final double EPSILON = 1e-4;

    @Test
    public void testFileLoading() throws IOException, EmptyFileException {
        final SortedSet<HmfGenomeRegion> regions = HmfGenomeFileLoader.fromFile(BASE_PATH + "test_gene_panel.tsv").get("13");
        assertEquals(2, regions.size());
        assertEquals(27, regions.first().exome().size());
        assertEquals(7, regions.last().exome().size());
    }

}
