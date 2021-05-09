package com.hartwig.hmftools.lilac.misc;

import static com.hartwig.hmftools.lilac.LilacConstants.HLA_TRANSCRIPTS;
import static com.hartwig.hmftools.lilac.LilacConstants.LOCI_POSITION;
import static com.hartwig.hmftools.lilac.ReferenceData.populateHlaTranscripts;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.lilac.LociPosition;

import org.junit.Test;

public class LociPositionTest
{
//    private val transcripts = HmfGenePanelSupplier.allGenesMap37()
//    private val aTranscript = transcripts[LilacApplication.HLA_A]!!
//    private val bTranscript = transcripts[LilacApplication.HLA_B]!!
//    private val hlaTranscripts = listOf(aTranscript, bTranscript, transcripts[LilacApplication.HLA_C]!!)

    @Test
    public void testA()
    {
        populateHlaTranscripts();

        HmfTranscriptRegion aTranscript = HLA_TRANSCRIPTS.get(0);
        int minLoci = LOCI_POSITION.nucelotideLoci((int)aTranscript.codingStart());
        assertTrue(minLoci >= 0);

        int maxLoci = LOCI_POSITION.nucelotideLoci((int)aTranscript.codingEnd());

        for (int i = minLoci; i < maxLoci; ++i)
        {
            int position = LOCI_POSITION.forwardPosition(i, aTranscript);
            int loci = LOCI_POSITION.nucelotideLoci(position);
            assertEquals(i, loci);
        }
    }

    @Test
    public void testB()
    {
        populateHlaTranscripts();
        HmfTranscriptRegion bTranscript = HLA_TRANSCRIPTS.get(1);

        int minLoci = LOCI_POSITION.nucelotideLoci((int)bTranscript.codingEnd());
        int maxLoci = LOCI_POSITION.nucelotideLoci((int)bTranscript.codingStart());
        assertTrue(minLoci >= 0);
        assertTrue(maxLoci > minLoci);

        for (int i = minLoci; i < maxLoci; ++i)
        {
            int position = LOCI_POSITION.reversePosition(i, bTranscript);
            int loci = LOCI_POSITION.nucelotideLoci(position);
            assertEquals(i, loci);
        }
    }

}
