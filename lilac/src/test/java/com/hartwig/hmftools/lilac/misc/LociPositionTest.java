package com.hartwig.hmftools.lilac.misc;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_A;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_B;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_C;
import static com.hartwig.hmftools.lilac.LilacUtils.calcNucelotideLocus;
import static com.hartwig.hmftools.lilac.ReferenceData.populateHlaTranscripts;

import static junit.framework.TestCase.assertEquals;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.lilac.MhcClass;

import org.junit.Test;

public class LociPositionTest
{
    @Test
    public void testNucleotideLocus()
    {
        Map<String,TranscriptData> hlaTranscriptMap = Maps.newHashMap();
        populateHlaTranscripts(hlaTranscriptMap, V37, MhcClass.CLASS_1);

        List<TranscriptData> transcripts = hlaTranscriptMap.values().stream().toList();

        TranscriptData aTranscript = hlaTranscriptMap.get(HLA_A);

        int locus = calcNucelotideLocus(transcripts, aTranscript.CodingStart);
        assertEquals(0, locus);

        locus = calcNucelotideLocus(transcripts, aTranscript.CodingEnd);
        assertEquals(1097, locus);

        TranscriptData bTranscript = hlaTranscriptMap.get(HLA_B);

        locus = calcNucelotideLocus(transcripts, bTranscript.CodingEnd);
        assertEquals(0, locus);

        locus = calcNucelotideLocus(transcripts, bTranscript.CodingStart);
        assertEquals(1088, locus);

        TranscriptData cTranscript = hlaTranscriptMap.get(HLA_C);

        locus = calcNucelotideLocus(transcripts, cTranscript.CodingEnd);
        assertEquals(0, locus);

        locus = calcNucelotideLocus(transcripts, cTranscript.CodingStart);
        assertEquals(1100, locus);

        hlaTranscriptMap.clear();
        populateHlaTranscripts(hlaTranscriptMap, V38, MhcClass.CLASS_1);

        transcripts = hlaTranscriptMap.values().stream().toList();

        aTranscript = hlaTranscriptMap.get(HLA_A);

        locus = calcNucelotideLocus(transcripts, aTranscript.CodingStart);
        assertEquals(0, locus);

        locus = calcNucelotideLocus(transcripts, aTranscript.CodingEnd);
        assertEquals(1097, locus);

        bTranscript = hlaTranscriptMap.get(HLA_B);

        locus = calcNucelotideLocus(transcripts, bTranscript.CodingEnd);
        assertEquals(0, locus);

        locus = calcNucelotideLocus(transcripts, bTranscript.CodingStart);
        assertEquals(1088, locus);

        cTranscript = hlaTranscriptMap.get(HLA_C);

        locus = calcNucelotideLocus(transcripts, cTranscript.CodingEnd);
        assertEquals(0, locus);

        locus = calcNucelotideLocus(transcripts, cTranscript.CodingStart);
        assertEquals(1100, locus);

    }
}
