package com.hartwig.hmftools.common.genepanel;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.common.gene.TranscriptRegion;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;

import org.junit.Test;

public class HmfGenePanelSupplierTest {

    @Test
    public void canLoadGeneRegionsFromFile() throws IOException {
        final Set<String> panel = HmfGenePanelSupplier.hmfPanelGeneSet();
        final List<HmfGenomeRegion> geneRegions = HmfGenePanelSupplier.hmfPanelGeneList();
        assertEquals(panel.size(), geneRegions.size());
    }

    @Test
    public void loadedRegionsAreSortedCorrectly() {
        final SortedSetMultimap<String, HmfGenomeRegion> geneRegions = HmfGenePanelSupplier.allGeneMap();
        for (final String chromosome : geneRegions.keySet()) {
            long start = 0;
            for (final HmfGenomeRegion hmfGenomeRegion : geneRegions.get(chromosome)) {
                assertTrue(hmfGenomeRegion.start() >= start);
                start = hmfGenomeRegion.start();
            }
        }
    }

    @Test
    public void testManuallyAddedGenes()  {
        final List<String> allGenes = HmfGenePanelSupplier.allGeneList().stream().map(TranscriptRegion::gene).collect(Collectors.toList());
        assertTrue(allGenes.contains("C11orf95"));
        assertTrue(allGenes.contains("CDKN2Ap14ARF"));
    }

}
