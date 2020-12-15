package com.hartwig.hmftools.common.fusion;

import static java.lang.Math.abs;

import java.util.List;

import com.hartwig.hmftools.common.ensemblcache.ExonData;

public class TranscriptUtils
{
    public static final int CODING_BASES = 0;
    public static final int TOTAL_CODING_BASES = 1;

    public static int codingBasesToPhase(int codingBases)
    {
        // subtract 1 to get back to phasing starting at zero, ie start with the first coding base where position == coding start:
        // coding base: 1 2 3 4 5 6 7 8 9 10
        // phase:       0 1 2 0 1 2 0 1 2 0
        codingBases -= 1;

        if(codingBases >= 0)
            return codingBases % 3;

        // coding base: -6 -5 -4 -3 -2 -1
        // phase:        0  1  2  0  1  2

        int mod = abs(codingBases) % 3;

        if(mod == 0)
            return 0;

        return mod == 1 ? 2 : 1;
    }

    public static int tickPhaseForward(int phase) { return (phase + 1) % 3; }

    public static int calcPositionPhasing(final Transcript transcript, boolean isUpstream)
    {
        int codingBases = transcript.codingBases();

        // factor in insert sequence for the upstream partner
        if(isUpstream && !transcript.gene().insertSequence().isEmpty())
        {
            codingBases += transcript.gene().insertSequence().length();
        }

        return codingBasesToPhase(codingBases);
    }

    public static int[] calcCodingBases(int codingStart, int codingEnd, final List<ExonData> exonList, int position)
    {
        boolean inCodingRegion = false;
        boolean codingRegionEnded = false;

        int[] codingData = {0, 0};

        for (ExonData exonData : exonList)
        {
            int exonStart = exonData.ExonStart;
            int exonEnd = exonData.ExonEnd;

            if (!inCodingRegion)
            {
                if (exonEnd >= codingStart)
                {
                    // coding region begins in this exon
                    inCodingRegion = true;

                    codingData[TOTAL_CODING_BASES] += exonEnd - codingStart + 1;

                    // check whether the position falls in this exon and if so before or after the coding start
                    if (position >= codingStart)
                    {
                        if (position < exonEnd)
                            codingData[CODING_BASES] += position - codingStart + 1;
                        else
                            codingData[CODING_BASES] += exonEnd - codingStart + 1;
                    }
                }
            }
            else if (!codingRegionEnded)
            {
                if (exonStart > codingEnd)
                {
                    codingRegionEnded = true;
                }
                else if (exonEnd > codingEnd)
                {
                    // coding region ends in this exon
                    codingRegionEnded = true;

                    codingData[TOTAL_CODING_BASES] += codingEnd - exonStart + 1;

                    if (position >= exonStart)
                    {
                        if (position < codingEnd)
                            codingData[CODING_BASES] += position - exonStart + 1;
                        else
                            codingData[CODING_BASES] += codingEnd - exonStart + 1;
                    }
                }
                else
                {
                    // take all of the exon's bases
                    codingData[TOTAL_CODING_BASES] += exonEnd - exonStart + 1;

                    if (position >= exonStart)
                    {
                        if (position < exonEnd)
                            codingData[CODING_BASES] += position - exonStart + 1;
                        else
                            codingData[CODING_BASES] += exonEnd - exonStart + 1;
                    }
                }
            }
        }

        return codingData;
    }

}
