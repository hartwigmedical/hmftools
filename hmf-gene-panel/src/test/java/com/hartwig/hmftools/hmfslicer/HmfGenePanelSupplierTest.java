package com.hartwig.hmftools.hmfslicer;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;

import org.junit.Test;

public class HmfGenePanelSupplierTest {

    @Test
    public void canLoadGeneRegionsFromFile() throws IOException, EmptyFileException {
        final List<String> genes = HmfGenePanelBuilder.readGeneList();
        final List<HmfGenomeRegion> geneRegions = HmfGenePanelSupplier.asList();
        assertEquals(genes.size(), geneRegions.size());
    }

    @Test
    public void loadedRegionsAreSortedCorrectly() throws IOException, EmptyFileException {
        final SortedSetMultimap<String, HmfGenomeRegion> geneRegions = HmfGenePanelSupplier.asMap();
        for (final String chromosome : geneRegions.keySet()) {
            long start = 0;
            for (final HmfGenomeRegion hmfGenomeRegion : geneRegions.get(chromosome)) {
                assertTrue(hmfGenomeRegion.start() >= start);
                start = hmfGenomeRegion.start();
            }
        }
    }
}
