package com.hartwig.hmftools.pavereverse.dna;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.ReversePaveTestBase;
import com.hartwig.hmftools.pavereverse.gene.GeneTranscript;

import org.junit.Test;

public class InIntronUpstreamOfCodingStartTest extends ReversePaveTestBase
{
    @Test
    public void locationTest()
    {
        int[] exonStarts = { 10, 30, 50, 70, 90 };
        int codingStart = 55;
        int codingEnd = 94;
        // 12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
        // 1        10        20        30        40        50        60        70        80        90        100       110
        // ____-----++++++++++----------++++++++++----------+++++*****----------**********----------*****+++++------_____

        GeneData geneData = createEnsemblGeneData("id132", "BLAH", "1",  POS_STRAND, 5, 105);
        TranscriptData transcript = createTransExons(geneData.GeneId, 123, POS_STRAND, exonStarts, 9, codingStart, codingEnd, false, "whatever");
        GeneTranscript gt = new GeneTranscript(geneData, transcript);

        assertEquals(49, new InIntronUpstreamOfCodingStart(-5, -1).toStrandLocation(gt));
        assertEquals(48, new InIntronUpstreamOfCodingStart(-5, -2).toStrandLocation(gt));
        assertEquals(40, new InIntronUpstreamOfCodingStart(-6, 1).toStrandLocation(gt));
        assertEquals(41, new InIntronUpstreamOfCodingStart(-6, 2).toStrandLocation(gt));
        assertEquals(29, new InIntronUpstreamOfCodingStart(-15, -1).toStrandLocation(gt));
        assertEquals(20, new InIntronUpstreamOfCodingStart(-16, 1).toStrandLocation(gt));
        assertEquals(9, new InIntronUpstreamOfCodingStart(-25, -1).toStrandLocation(gt));
    }
}
