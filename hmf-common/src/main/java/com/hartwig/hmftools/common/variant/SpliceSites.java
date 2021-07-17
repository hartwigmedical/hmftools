package com.hartwig.hmftools.common.variant;

import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.Strand;

public class SpliceSites
{
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

    private static boolean isDonorWithOffset(final HmfTranscriptRegion transcript, int offset, long position)
    {
        // returns true if the position matches the donor site of the exon (where offset = 0 means first base after the exon)
        if(transcript.codingStart() == 0)
            return false;

        if(position < transcript.codingStart() || position > transcript.codingEnd())
            return false;

        if(transcript.strand() == Strand.FORWARD)
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

        if(transcript.strand() == Strand.FORWARD)
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
