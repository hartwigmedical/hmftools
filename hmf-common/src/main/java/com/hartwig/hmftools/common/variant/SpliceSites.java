package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.genome.region.Strand.FORWARD;

import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.Strand;

public final class SpliceSites
{
    public static int getDonorPosition(int position, int exonEnd, Strand strand)
    {
        // zero is the first base after the exon, so pos of 101, exon end of 100, is position zero
        if(strand == FORWARD)
            return position - exonEnd - 1;
        else
            return exonEnd - position - 1; // pos of 99, exon end of 100, is position zero
    }

    public static int getAcceptorPosition(int position, int exonEnd, Strand strand)
    {
        // zero is the first base before the exon, so pos of 99, exon end of 100, is position zero
        if(strand == FORWARD)
            return exonEnd - position - 1;
        else
            return position - exonEnd - 1; // pos of 101, exon end of 100, is position zero
    }

    public static boolean isDonorMinusOne(final HmfTranscriptRegion transcript, long position)
    {
        // check if the position matches the last exon base at a donor site
        return isDonorWithOffset(transcript, -1, position);
    }

    public static boolean isDonorPlusFive(final HmfTranscriptRegion transcript, long position)
    {
        // check if the position matches the 5th base into the intron after the exon at the donor site
        return isDonorWithOffset(transcript, 4, position);
    }

    public static boolean isAcceptorPlusThree(final HmfTranscriptRegion transcript, long position)
    {
        // check if the position matches the 3rd base into the intron before the exon at the acceptor site
        return isAcceptorWithOffset(transcript, 2, position);
    }

    // TODO - combine these methods when HmfTranscriptRegion is removed
    private static boolean isDonorWithOffset(final HmfTranscriptRegion transcript, int offset, long position)
    {
        // returns true if the position matches the donor site of the exon (where offset = 0 means first base after the exon)
        if(transcript.codingStart() == 0)
            return false;

        if(position < transcript.codingStart() || position > transcript.codingEnd())
            return false;

        if(transcript.strand() == FORWARD)
        {
            for(int i = 0; i < transcript.exome().size() - 1; i++)
            {
                long donorSite = transcript.exome().get(i).end() + 1;

                if(position == donorSite + offset)
                    return true;
            }
        }
        else
        {
            for(int i = 1; i < transcript.exome().size(); i++)
            {
                long donorSite = transcript.exome().get(i).start() - 1;

                if(position == donorSite - offset)
                    return true;
            }
        }

        return false;
    }

    private static boolean isAcceptorWithOffset(final HmfTranscriptRegion transcript, int offset, long position)
    {
        // returns true if the position matches the acceptor site of the exon (where offset = 0 means first base before the exon)
        if(transcript.codingStart() == 0)
            return false;

        if(position < transcript.codingStart() || position > transcript.codingEnd())
            return false;

        if(transcript.strand() == FORWARD)
        {
            for(int i = 1; i < transcript.exome().size(); i++)
            {
                long acceptorSide = transcript.exome().get(i).start() - 1;

                if(position == acceptorSide - offset)
                    return true;
            }
        }
        else
        {
            for(int i = 0; i < transcript.exome().size() - 1; i++)
            {
                long donorSite = transcript.exome().get(i).end() + 1;

                if(position == donorSite + offset)
                    return true;
            }
        }

        return false;
    }

}
