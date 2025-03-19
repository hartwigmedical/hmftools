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

public class InIntronDownstreamOfCodingEndTest extends ReversePaveTestBase
{
    @Test
    public void locationTest()
    {
        int[] exonStarts = { 10, 30, 50, 70, 90 };
        int codingStart = 35;
        int codingEnd = 74;
        // 12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
        // 1        10        20        30        40        50        60        70        80        90        100       110
        // ____-----++++++++++----------+++++*****----------**********----------*****+++++----------++++++++++------_____

        GeneData geneData = createEnsemblGeneData("id132", "BLAH", "1",  POS_STRAND, 5, 105);
        TranscriptData transcript = createTransExons(geneData.GeneId, 123, POS_STRAND, exonStarts, 9, codingStart, codingEnd, false, "whatever");
        GeneTranscript gt = new GeneTranscript(geneData, transcript);

        assertEquals(80, new InIntronDownstreamOfCodingEnd(5, 1).toStrandLocation(gt));
        assertEquals(81, new InIntronDownstreamOfCodingEnd(5, 2).toStrandLocation(gt));
        assertEquals(89, new InIntronDownstreamOfCodingEnd(6, -1).toStrandLocation(gt));
        assertEquals(88, new InIntronDownstreamOfCodingEnd(6, -2).toStrandLocation(gt));
        assertEquals(100, new InIntronDownstreamOfCodingEnd(15, 1).toStrandLocation(gt));
    }
}
