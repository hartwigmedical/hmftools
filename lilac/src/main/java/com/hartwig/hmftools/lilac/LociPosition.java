package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.common.gene.TranscriptUtils.calcCodingBases;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.bed.ImmutableNamedBed;
import com.hartwig.hmftools.common.genome.bed.NamedBed;

import java.util.List;

public class LociPosition
{
    private final List<TranscriptData> mTranscripts;

    public LociPosition(final List<TranscriptData> transcripts)
    {
        mTranscripts = transcripts;
    }

    public final int calcNucelotideLocus(int position)
    {
        for(TranscriptData transData : mTranscripts)
        {
            if(position < transData.CodingStart || position > transData.CodingEnd)
                continue;

            // locus is a zero-based index, so the first coding base has locus of 0
            return calcCodingBases(transData, position).CodingBases - 1;
        }
        return -1;
    }

    public static List<NamedBed> codingRegions(final String geneName, final String chromosome, final TranscriptData transcript)
    {
        final List<NamedBed> result = Lists.newArrayList();

        int codingStart = transcript.CodingStart;
        int codingEnd = transcript.CodingEnd;

        for (ExonData exon : transcript.exons())
        {
            if (codingStart <= exon.End && codingEnd >= exon.Start)
            {
                result.add(ImmutableNamedBed.builder()
                        .chromosome(chromosome)
                        .start(Math.max(codingStart, exon.Start))
                        .end(Math.min(codingEnd, exon.End))
                        .name(geneName)
                        .build());
            }
        }

        return result;
    }
}
