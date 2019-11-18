package com.hartwig.hmftools.common.genome.genepanel;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.BEDFileLoader;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;

import org.junit.Ignore;
import org.junit.Test;

public class HmfExonPanelBedTest {

    @Test
    public void testExonPlusMinusBases() {
        HmfTranscriptRegion firstTranscript = HmfGenePanelSupplier.allGeneList37().get(0);
        HmfExonRegion firstExon = firstTranscript.exome().get(0);
        List<GenomeRegion> regions = HmfExonPanelBed.createRegions(Lists.newArrayList(firstTranscript, firstTranscript));
        assertEquals(firstTranscript.exome().size(), regions.size());

        GenomeRegion region = regions.get(0);
        assertEquals(firstExon.start() - HmfExonPanelBed.EXTRA_BASES, region.start());
        assertEquals(firstExon.end() + HmfExonPanelBed.EXTRA_BASES, region.end());
    }

    @Test
    @Ignore
    public void write19() throws IOException {
        String filename = "/Users/jon/hmf/resources/ActionableExonPanel.hg19.bed";
        HmfExonPanelBed.write19File(filename);
        BEDFileLoader.fromBedFile(filename);
    }

    @Test
    @Ignore
    public void write38() throws IOException {
        String filename = "/Users/jon/hmf/resources/ActionableExonPanel.hg38.bed";
        HmfExonPanelBed.write38File(filename);
        BEDFileLoader.fromBedFile(filename);
    }
}
