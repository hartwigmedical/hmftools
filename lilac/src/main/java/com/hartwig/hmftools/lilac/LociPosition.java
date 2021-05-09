package com.hartwig.hmftools.lilac;

import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.Strand;

import java.util.List;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public class LociPosition
{
    private final List<HmfTranscriptRegion> mTranscripts;

    public LociPosition()
    {
        mTranscripts = Lists.newArrayList();
    }

    public void initialise(final List<HmfTranscriptRegion> transcripts)
    {
        mTranscripts.addAll(transcripts);
    }

    public final int nucelotideLoci(int position)
    {
        for(HmfTranscriptRegion transcript : mTranscripts)
        {
            if((long) position < transcript.start() || (long) position > transcript.end())
            {
                continue;
            }
            if(transcript.strand() == Strand.FORWARD)
            {
                return forwardLoci(position, transcript);
            }
            return reverseLoci(position, transcript);
        }
        return -1;
    }

    /* UNUSED - confirm then delete
    public final List<Integer> position(int codingLoci)
    {
        List result = new ArrayList();
        for(HmfTranscriptRegion transcript : mTranscripts)
        {
            if(transcript.strand() == Strand.FORWARD)
            {
                result.add(forwardPosition(codingLoci, transcript));
                continue;
            }
            result.add(reversePosition(codingLoci, transcript));
        }
        return result;
    }

    public final int position(int codingLoci, @NotNull HmfTranscriptRegion transcript)
    {
        return transcript.strand() == Strand.FORWARD
                ? forwardPosition(codingLoci, transcript)
                : reversePosition(codingLoci, transcript);
    }
     */

    public final int reverseLoci(int position, final HmfTranscriptRegion transcript)
    {
        if(transcript.strand() != Strand.REVERSE)
            return -1;

        long codingStart = transcript.codingStart();
        long codingEnd = transcript.codingEnd();
        int currentLoci = 0;

        for(int i = transcript.exome().size() - 1; i >= 0; i--)
        {
            HmfExonRegion exon = transcript.exome().get(i);

            if(exon.end() < codingStart || exon.start() > codingEnd)
                continue;

            long exonStartPosition = Math.max(codingStart, exon.start());
            long exonEndPosition = Math.min(codingEnd, exon.end());
            int exonLength = (int) (exonEndPosition - exonStartPosition + 1L);
            int exonStartLoci = currentLoci;
            int exonEndLoci = exonStartLoci + exonLength - 1;
            long l = position;
            if(exonStartPosition <= l && exonEndPosition >= l)
            {
                return (int) ((long) (currentLoci - position) + exonEndPosition);
            }
            currentLoci = exonEndLoci + 1;
        }
        return -1;
    }

    public final int forwardLoci(int position, final HmfTranscriptRegion transcript)
    {
        if(transcript.strand() != Strand.FORWARD)
            return -1;

        long codingStart = transcript.codingStart();
        long codingEnd = transcript.codingEnd();
        int currentLoci = 0;
        for(HmfExonRegion exon : transcript.exome())
        {
            if(exon.end() < codingStart || exon.start() > codingEnd)
                continue;

            long exonStartPosition = Math.max(codingStart, exon.start());
            long exonEndPosition = Math.min(codingEnd, exon.end());
            int exonLength = (int) (exonEndPosition - exonStartPosition + 1L);
            int exonStartLoci = currentLoci;
            int exonEndLoci = exonStartLoci + exonLength - 1;
            long l = position;
            if(exonStartPosition <= l && exonEndPosition >= l)
            {
                return (int) ((long) (currentLoci + position) - exonStartPosition);
            }
            currentLoci = exonEndLoci + 1;
        }
        return -1;
    }

    public final int reversePosition(int codingLoci, final HmfTranscriptRegion transcript)
    {
        if(transcript.strand() != Strand.REVERSE)
            return -1;

        long codingStart = transcript.codingStart();
        long codingEnd = transcript.codingEnd();
        int currentLoci = 0;

        for(int i = transcript.exome().size() - 1; i >= 0; i--)
        {
            HmfExonRegion exon = transcript.exome().get(i);

            int exonEndLoci = 0;
            if(exon.end() < codingStart || exon.start() > codingEnd)
                continue;

            long exonStartPosition = Math.max(codingStart, exon.start());
            long exonEndPosition = Math.min(codingEnd, exon.end());
            int exonLength = (int) (exonEndPosition - exonStartPosition + 1L);
            int exonStartLoci = currentLoci;
            int n = codingLoci;
            if(exonStartLoci <= n && (exonEndLoci = exonStartLoci + exonLength - 1) >= n)
            {
                return (int) (exonEndPosition - (long) codingLoci + (long) exonStartLoci);
            }
            currentLoci = exonEndLoci + 1;
        }
        return -1;
    }

    public final int forwardPosition(int codingLoci, HmfTranscriptRegion transcript)
    {
        if(transcript.strand() != Strand.FORWARD)
            return -1;

        long codingStart = transcript.codingStart();
        long codingEnd = transcript.codingEnd();
        int currentLoci = 0;
        for(HmfExonRegion exon : transcript.exome())
        {
            int exonEndLoci = 0;

            if(exon.end() < codingStart || exon.start() > codingEnd)
                continue;

            long exonStartPosition = Math.max(codingStart, exon.start());
            long exonEndPosition = Math.min(codingEnd, exon.end());
            int exonLength = (int) (exonEndPosition - exonStartPosition + 1L);
            int exonStartLoci = currentLoci;
            int n = codingLoci;
            if(exonStartLoci <= n && (exonEndLoci = exonStartLoci + exonLength - 1) >= n)
            {
                return (int) (exonStartPosition + (long) codingLoci - (long) exonStartLoci);
            }
            currentLoci = exonEndLoci + 1;
        }
        return -1;
    }
}
