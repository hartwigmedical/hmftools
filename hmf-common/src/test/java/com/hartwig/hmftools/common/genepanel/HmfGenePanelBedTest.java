package com.hartwig.hmftools.common.genepanel;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;

import org.junit.Ignore;
import org.junit.Test;

public class HmfGenePanelBedTest {

    @Test
    public void testExonPlusMinusBases() {
        HmfTranscriptRegion firstRegion = HmfGenePanelSupplier.allGeneList37().get(0);
        List<GenomeRegion> regions = HmfGenePanelBed.createRegions(Lists.newArrayList(firstRegion, firstRegion));
        assertEquals(1, regions.size());

        GenomeRegion region = regions.get(0);
        assertEquals(firstRegion.start() - HmfGenePanelBed.EXTRA_BASES, region.start());
        assertEquals(firstRegion.end() + HmfGenePanelBed.EXTRA_BASES, region.end());
    }


    @Test
    @Ignore
    public void write19() throws IOException {
        HmfGenePanelBed.write37File("/Users/jon/hmf/resources/GenePanel.hg19.bed");
    }

    @Ignore
    public void write38() throws IOException {
        HmfGenePanelBed.write38File("~/hmf/resources/GenePanel.hg38.bed");
    }


}
