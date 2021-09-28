package com.hartwig.hmftools.common.genome.genepanel;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class HmfExonPanelBedTest {

    @Test
    public void testReverseTranscript() {
        final HmfTranscriptRegion transcript = HmfGenePanelSupplier.allGenesMap37().get("TP53");
        HmfExonRegion firstExon = transcript.exons().get(0);
        HmfExonRegion secondExon = transcript.exons().get(1);
        HmfExonRegion thirdExon = transcript.exons().get(2);
        HmfExonRegion finalCodingExon = transcript.exons().get(transcript.exons().size() - 2);

        List<? extends GenomeRegion> regions =
                HmfExonPanelBed.createNamedCodingRegions(false, Sets.newHashSet("TP53"), Lists.newArrayList(transcript, transcript));

        assertRegion(regions.get(0), transcript.codingStart(), firstExon.end() + 10);
        assertRegion(regions.get(1), secondExon.start() - 10, secondExon.end() + 10);
        assertRegion(regions.get(2), thirdExon.start() - 10, thirdExon.end() + 10);

        assertRegion(regions.get(9), finalCodingExon.start() - 10, transcript.codingEnd());
    }

    private void assertRegion(@NotNull final GenomeRegion victim, long expectedStart, long expectedEnd) {
        assertEquals(expectedStart, victim.start());
        assertEquals(expectedEnd, victim.end());

    }

}
