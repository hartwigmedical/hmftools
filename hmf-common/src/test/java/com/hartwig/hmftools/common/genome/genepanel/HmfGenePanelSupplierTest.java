package com.hartwig.hmftools.common.genome.genepanel;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.TranscriptRegion;

import org.junit.Test;

public class HmfGenePanelSupplierTest {

    @Test
    public void allRegionsAreSortedCorrectly() {
        final SortedSetMultimap<String, HmfTranscriptRegion> geneRegions = HmfGenePanelSupplier.allGenesPerChromosomeMap37();
        for (final String chromosome : geneRegions.keySet()) {
            int start = 0;
            for (final HmfTranscriptRegion hmfTranscriptRegion : geneRegions.get(chromosome)) {
                assertTrue(hmfTranscriptRegion.start() >= start);
                start = hmfTranscriptRegion.start();
            }
        }
    }

    @Test
    public void allGenesAreUnique37() {
        List<HmfTranscriptRegion> regionList = HmfGenePanelSupplier.allGeneList37();
        Map<String, HmfTranscriptRegion> regionsByGene = HmfGenePanelSupplier.allGenesMap37();

        assertEquals(regionList.size(), regionsByGene.values().size());
    }

    @Test
    public void allGenesAreUnique38() {
        List<HmfTranscriptRegion> regionList = HmfGenePanelSupplier.allGeneList38();
        Map<String, HmfTranscriptRegion> regionsByGene = HmfGenePanelSupplier.allGenesMap38();

        assertEquals(regionList.size(), regionsByGene.values().size());
    }

    @Test
    public void verifyManuallyAddedGenesArePresent() {
        final List<String> allGenes =
                HmfGenePanelSupplier.allGeneList37().stream().map(TranscriptRegion::geneName).collect(Collectors.toList());
        assertTrue(allGenes.contains("C11orf95"));
        assertTrue(allGenes.contains("CDKN2Ap14ARF"));
    }
}
