package com.hartwig.hmftools.lilac.misc;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_A;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_B;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_C;
import static com.hartwig.hmftools.lilac.ReferenceData.populateHlaTranscripts;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.lilac.LociPosition;

import org.junit.Test;

public class LociPositionTest
{
    @Test
    public void testNucleotideLocus()
    {
        Map<String, TranscriptData> hlaTranscriptMap = Maps.newHashMap();
        populateHlaTranscripts(hlaTranscriptMap, V37);
        LociPosition lociPosition = new LociPosition(hlaTranscriptMap.values().stream().collect(Collectors.toList()));

        TranscriptData aTranscript = hlaTranscriptMap.get(HLA_A);

        int locus = lociPosition.calcNucelotideLocus(aTranscript.CodingStart);
        assertEquals(0, locus);

        locus = lociPosition.calcNucelotideLocus(aTranscript.CodingEnd);
        assertEquals(1097, locus);

        TranscriptData bTranscript = hlaTranscriptMap.get(HLA_B);

        locus = lociPosition.calcNucelotideLocus(bTranscript.CodingEnd);
        assertEquals(0, locus);

        locus = lociPosition.calcNucelotideLocus(bTranscript.CodingStart);
        assertEquals(1088, locus);

        TranscriptData cTranscript = hlaTranscriptMap.get(HLA_C);

        locus = lociPosition.calcNucelotideLocus(cTranscript.CodingEnd);
        assertEquals(0, locus);

        locus = lociPosition.calcNucelotideLocus(cTranscript.CodingStart);
        assertEquals(1100, locus);

        hlaTranscriptMap.clear();
        populateHlaTranscripts(hlaTranscriptMap, V38);
        lociPosition = new LociPosition(hlaTranscriptMap.values().stream().collect(Collectors.toList()));

        aTranscript = hlaTranscriptMap.get(HLA_A);

        locus = lociPosition.calcNucelotideLocus(aTranscript.CodingStart);
        assertEquals(0, locus);

        locus = lociPosition.calcNucelotideLocus(aTranscript.CodingEnd);
        assertEquals(1097, locus);

        bTranscript = hlaTranscriptMap.get(HLA_B);

        locus = lociPosition.calcNucelotideLocus(bTranscript.CodingEnd);
        assertEquals(0, locus);

        locus = lociPosition.calcNucelotideLocus(bTranscript.CodingStart);
        assertEquals(1088, locus);

        cTranscript = hlaTranscriptMap.get(HLA_C);

        locus = lociPosition.calcNucelotideLocus(cTranscript.CodingEnd);
        assertEquals(0, locus);

        locus = lociPosition.calcNucelotideLocus(cTranscript.CodingStart);
        assertEquals(1100, locus);

    }
}
