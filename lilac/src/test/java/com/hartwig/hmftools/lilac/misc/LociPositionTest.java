package com.hartwig.hmftools.lilac.misc;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_A;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_B;
import static com.hartwig.hmftools.lilac.ReferenceData.populateHlaTranscripts;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.lilac.LociPosition;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import org.junit.Test;

public class LociPositionTest
{
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
        Map<String, TranscriptData> hlaTranscriptMap = Maps.newHashMap();
        populateHlaTranscripts(hlaTranscriptMap, V37);
        LociPosition lociPosition = new LociPosition(hlaTranscriptMap.values().stream().collect(Collectors.toList()));

        TranscriptData aTranscript = hlaTranscriptMap.get(HLA_A);
        int minLoci = lociPosition.nucelotideLoci((int)aTranscript.CodingStart);
        assertTrue(minLoci >= 0);

        int maxLoci = lociPosition.nucelotideLoci((int)aTranscript.CodingEnd);

        for (int i = minLoci; i < maxLoci; ++i)
        {
            int position = lociPosition.forwardPosition(i, aTranscript);
            int loci = lociPosition.nucelotideLoci(position);
            assertEquals(i, loci);
        }
    }

    @Test
    public void testB()
    {
        Map<String,TranscriptData> hlaTranscriptMap = Maps.newHashMap();
        populateHlaTranscripts(hlaTranscriptMap, V37);
        TranscriptData bTranscript = hlaTranscriptMap.get(HLA_B);
        LociPosition lociPosition = new LociPosition(hlaTranscriptMap.values().stream().collect(Collectors.toList()));

        int minLoci = lociPosition.nucelotideLoci((int)bTranscript.CodingEnd);
        int maxLoci = lociPosition.nucelotideLoci((int)bTranscript.CodingStart);
        assertTrue(minLoci >= 0);
        assertTrue(maxLoci > minLoci);

        for (int i = minLoci; i < maxLoci; ++i)
        {
            int position = lociPosition.reversePosition(i, bTranscript);
            int loci = lociPosition.nucelotideLoci(position);
            assertEquals(i, loci);
        }
    }

}
