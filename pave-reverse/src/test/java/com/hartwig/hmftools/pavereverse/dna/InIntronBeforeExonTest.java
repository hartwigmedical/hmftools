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

public class InIntronBeforeExonTest extends ReversePaveTestBase
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

        assertEquals(30 - 1, new InIntronBeforeExon(7, 1).toStrandLocation(gt));
        assertEquals(30 - 2, new InIntronBeforeExon(7, 2).toStrandLocation(gt));
        assertEquals(50 - 1, new InIntronBeforeExon(18, 1).toStrandLocation(gt));
        assertEquals(50 - 3, new InIntronBeforeExon(18, 3).toStrandLocation(gt));
    }
}
