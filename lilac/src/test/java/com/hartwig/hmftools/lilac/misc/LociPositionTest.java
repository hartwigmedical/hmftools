package com.hartwig.hmftools.lilac.misc;

import static com.hartwig.hmftools.lilac.LilacConstants.HLA_TRANSCRIPTS;
import static com.hartwig.hmftools.lilac.LilacConstants.LOCI_POSITION;
import static com.hartwig.hmftools.lilac.ReferenceData.populateHlaTranscripts;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.lilac.LociPosition;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import org.junit.Test;

public class LociPositionTest
{
//    private val transcripts = HmfGenePanelSupplier.allGenesMap37()
//    private val aTranscript = transcripts[LilacApplication.HLA_A]!!
//    private val bTranscript = transcripts[LilacApplication.HLA_B]!!
//    private val hlaTranscripts = listOf(aTranscript, bTranscript, transcripts[LilacApplication.HLA_C]!!)

    @Test
    public void testAlleleMaps()
    {
        Map<HlaAllele,Integer> map = Maps.newHashMap();

        HlaAllele allele1 = HlaAllele.fromString("A*01:01");
        HlaAllele allele2 = HlaAllele.fromString("A*01:01");
        HlaAllele allele3 = HlaAllele.fromString("A*01:01:01");
        HlaAllele allele4 = HlaAllele.fromString("A*01:02");

        map.put(allele1, 1);
        map.put(allele2, 2);
        map.put(allele3, 3);
        map.put(allele4, 4);

        assertEquals(map.size(), 3);

        Integer val1 = map.get(allele1);
        Integer val2 = map.get(allele2);

        assertEquals(val1, val2);

        map.remove(allele1);
        assertEquals(map.size(), 2);

        map.remove(allele3);
        assertEquals(map.size(), 1);

        map.remove(allele4);
        assertTrue(map.isEmpty());
    }

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
