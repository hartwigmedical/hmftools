package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
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

    public final int nucelotideLoci(int position)
    {
        for(TranscriptData transcript : mTranscripts)
        {
            if(position < transcript.TransStart || position > transcript.TransEnd)
            {
                continue;
            }
            if(transcript.Strand == POS_STRAND)
            {
                return forwardLoci(position, transcript);
            }
            return reverseLoci(position, transcript);
        }
        return -1;
    }

    public final int reverseLoci(int position, final TranscriptData transcript)
    {
        if(transcript.Strand != NEG_STRAND)
            return -1;

        int codingStart = transcript.CodingStart;
        int codingEnd = transcript.CodingEnd;
        int currentLoci = 0;

        for(int i = transcript.exons().size() - 1; i >= 0; i--)
        {
            ExonData exon = transcript.exons().get(i);

            if(exon.End < codingStart || exon.Start > codingEnd)
                continue;

            int exonStartPosition = Math.max(codingStart, exon.Start);
            int exonEndPosition = Math.min(codingEnd, exon.End);
            int exonLength = exonEndPosition - exonStartPosition + 1;
            int exonStartLoci = currentLoci;
            int exonEndLoci = exonStartLoci + exonLength - 1;
            int l = position;
            if(exonStartPosition <= l && exonEndPosition >= l)
            {
                return currentLoci - position + exonEndPosition;
            }
            currentLoci = exonEndLoci + 1;
        }
        return -1;
    }

    public final int forwardLoci(int position, final TranscriptData transcript)
    {
        if(transcript.Strand != POS_STRAND)
            return -1;

        int codingStart = transcript.CodingStart;
        int codingEnd = transcript.CodingEnd;
        int currentLoci = 0;
        for(ExonData exon : transcript.exons())
        {
            if(exon.End < codingStart || exon.Start > codingEnd)
                continue;

            int exonStartPosition = Math.max(codingStart, exon.Start);
            int exonEndPosition = Math.min(codingEnd, exon.End);
            int exonLength = exonEndPosition - exonStartPosition + 1;
            int exonStartLoci = currentLoci;
            int exonEndLoci = exonStartLoci + exonLength - 1;
            int l = position;
            if(exonStartPosition <= l && exonEndPosition >= l)
            {
                return currentLoci + position - exonStartPosition;
            }
            currentLoci = exonEndLoci + 1;
        }
        return -1;
    }

    public final int reversePosition(int codingLoci, final TranscriptData transcript)
    {
        if(transcript.Strand != NEG_STRAND)
            return -1;

        int codingStart = transcript.CodingStart;
        int codingEnd = transcript.CodingEnd;
        int currentLoci = 0;

        for(int i = transcript.exons().size() - 1; i >= 0; i--)
        {
            ExonData exon = transcript.exons().get(i);

            int exonEndLoci = 0;
            if(exon.End < codingStart || exon.Start > codingEnd)
                continue;

            int exonStartPosition = Math.max(codingStart, exon.Start);
            int exonEndPosition = Math.min(codingEnd, exon.End);
            int exonLength = exonEndPosition - exonStartPosition + 1;
            int exonStartLoci = currentLoci;
            int n = codingLoci;
            if(exonStartLoci <= n && (exonEndLoci = exonStartLoci + exonLength - 1) >= n)
            {
                return exonEndPosition - codingLoci + exonStartLoci;
            }
            currentLoci = exonEndLoci + 1;
        }
        return -1;
    }

    public final int forwardPosition(int codingLoci, TranscriptData transcript)
    {
        if(transcript.Strand != POS_STRAND)
            return -1;

        int codingStart = transcript.CodingStart;
        int codingEnd = transcript.CodingEnd;
        int currentLoci = 0;
        for(ExonData exon : transcript.exons())
        {
            int exonEndLoci = 0;

            if(exon.End < codingStart || exon.Start > codingEnd)
                continue;

            int exonStartPosition = Math.max(codingStart, exon.Start);
            int exonEndPosition = Math.min(codingEnd, exon.End);
            int exonLength = exonEndPosition - exonStartPosition + 1;
            int exonStartLoci = currentLoci;
            int n = codingLoci;
            if(exonStartLoci <= n && (exonEndLoci = exonStartLoci + exonLength - 1) >= n)
            {
                return exonStartPosition + codingLoci - exonStartLoci;
            }
            currentLoci = exonEndLoci + 1;
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
