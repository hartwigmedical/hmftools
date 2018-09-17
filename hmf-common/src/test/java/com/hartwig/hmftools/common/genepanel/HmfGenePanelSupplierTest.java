package com.hartwig.hmftools.common.genepanel;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.region.TranscriptRegion;

import org.junit.Test;

public class HmfGenePanelSupplierTest {

    @Test
    public void canLoadGenePanelFromFile() throws IOException {
        final Set<String> panel = HmfGenePanelSupplier.hmfPanelGeneSet();
        final List<HmfTranscriptRegion> geneRegions = HmfGenePanelSupplier.hmfPanelGeneList();
        assertEquals(panel.size(), geneRegions.size());
    }

    @Test
    public void allRegionsAreSortedCorrectly() {
        final SortedSetMultimap<String, HmfTranscriptRegion> geneRegions = HmfGenePanelSupplier.allGenesPerChromosomeMap();
        for (final String chromosome : geneRegions.keySet()) {
            long start = 0;
            for (final HmfTranscriptRegion hmfTranscriptRegion : geneRegions.get(chromosome)) {
                assertTrue(hmfTranscriptRegion.start() >= start);
                start = hmfTranscriptRegion.start();
            }
        }
    }

    @Test
    public void allGenesAreUnique() {
        List<HmfTranscriptRegion> regionList = HmfGenePanelSupplier.allGeneList();
        Map<String, HmfTranscriptRegion> regionsByGene = HmfGenePanelSupplier.allGenesMap();

        assertEquals(regionList.size(), regionsByGene.values().size());
    }

    @Test
    public void verifyManuallyAddedGenesArePresent()  {
        final List<String> allGenes = HmfGenePanelSupplier.allGeneList().stream().map(TranscriptRegion::gene).collect(Collectors.toList());
        assertTrue(allGenes.contains("C11orf95"));
        assertTrue(allGenes.contains("CDKN2Ap14ARF"));
    }
}
