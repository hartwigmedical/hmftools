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

public class InExonTest extends ReversePaveTestBase
{
    @Test
    public void locationTest()
    {
        int[] exonStarts = { 10, 30, 50 };
        int codingStart = 15;
        int codingEnd = 55;

        GeneData geneData = createEnsemblGeneData("id132", "BLAH", "1",  POS_STRAND, 5, 105);
        TranscriptData transcript = createTransExons(geneData.GeneId, 123, POS_STRAND, exonStarts, 10, codingStart, codingEnd, false, "whatever");
        GeneTranscript gt = new GeneTranscript(geneData, transcript);

        assertEquals(15, new InExon(1).toStrandLocation(gt));
        assertEquals(20, new InExon(6).toStrandLocation(gt));
        assertEquals(30, new InExon(7).toStrandLocation(gt));
        assertEquals(30 + 10, new InExon(17).toStrandLocation(gt));
        assertEquals(50, new InExon(18).toStrandLocation(gt));
        assertEquals(55, new InExon(23).toStrandLocation(gt));
    }
}
