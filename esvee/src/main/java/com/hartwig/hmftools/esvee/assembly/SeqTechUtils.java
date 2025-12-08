package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;

import static com.hartwig.hmftools.esvee.common.SvConstants.isUltima;

import java.util.List;

import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.read.ReadAdjustments;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;

public final class SeqTechUtils
{
    public static final int ULTIMA_READ_MISMATCH_MEDIUM_REPEAT_COUNT = 7;
    public static final int ULTIMA_READ_MISMATCH_LONG_REPEAT_COUNT = 11;

    public static final int SBX_PRIME_POSITION_RANGE_THRESHOLD = 6;

    public static void setSeqTechSpecifics()
    {
        if(isUltima())
        {
            AssemblyConstants.READ_MISMATCH_MEDIUM_REPEAT_COUNT = ULTIMA_READ_MISMATCH_MEDIUM_REPEAT_COUNT;
            AssemblyConstants.READ_MISMATCH_LONG_REPEAT_COUNT = ULTIMA_READ_MISMATCH_LONG_REPEAT_COUNT;
        }
    }

    public static void trimIlluminaAdapterBases(final Read read)
    {
        // trim any 3' bases extending past the unclipped 5' position
        if(!read.isPairedRead() || read.mateRead() == null)
            return;

        Read mateRead = read.mateRead();

        if(read.orientation() == mateRead.orientation())
            return;

        if(!read.chromosome().equals(mateRead.chromosome()))
            return;

        if(read.orientation().isForward())
        {
            int threePrimeEnd = read.unclippedEnd();
            int mateFivePrimeEnd = mateRead.unclippedEnd();

            if(threePrimeEnd > mateFivePrimeEnd)
            {
                read.trimBases(threePrimeEnd - mateFivePrimeEnd, false);
            }
        }
        else
        {
            int threePrimeStart = read.unclippedStart();
            int mateFivePrimeStart = mateRead.unclippedStart();

            if(threePrimeStart < mateFivePrimeStart)
            {
                read.trimBases(mateFivePrimeStart - threePrimeStart, true);
            }
        }
    }

    public static boolean trimSbxUncertainBases(final Read read)
    {
        // trim uncertain bases from both ends
        int lastIndex = read.basesLength() - 1;
        int trimCountStart = ReadAdjustments.findLowBaseQualTrimCount(read, 0, lastIndex, true, true);
        int trimCountEnd = ReadAdjustments.findLowBaseQualTrimCount(read, 0, lastIndex, false, true);

        if(trimCountStart > 0 && trimCountEnd > 0 && trimCountStart + trimCountEnd >= read.basesLength())
        {
            // apply the shorter of the 2 first then reassess
            if(trimCountStart < trimCountEnd)
            {
                read.trimBases(trimCountStart, true);

                lastIndex = read.basesLength() - 1;
                trimCountEnd = ReadAdjustments.findLowBaseQualTrimCount(read, 0, lastIndex, false, true);

                if(trimCountEnd > 0)
                    read.trimBases(trimCountEnd, false);
            }
            else
            {
                read.trimBases(trimCountEnd, false);

                lastIndex = read.basesLength() - 1;
                trimCountStart = ReadAdjustments.findLowBaseQualTrimCount(read, 0, lastIndex, false, true);

                if(trimCountStart > 0)
                    read.trimBases(trimCountStart, true);
            }
        }
        else
        {
            if(trimCountStart > 0)
                read.trimBases(trimCountStart, true);

            if(trimCountEnd > 0)
                read.trimBases(trimCountEnd, false);
        }

        return trimCountStart > 0 || trimCountEnd > 0;
    }

    public static boolean passSbxDistinctPrimePositionsFilter(final List<SupportRead> support)
    {
        // filter the 5' and 3' ranges are too tight
        int minFivePrimePos = -1;
        int maxFivePrimePos = -1;
        int minThreePrimePos = -1;
        int maxThreePrimePos = -1;

        for(SupportRead read : support)
        {
            int fivePrimePos, threePrimePos;

            if(read.orientation().isForward())
            {
                fivePrimePos = read.untrimmedStart();
                threePrimePos = read.untrimmedEnd();
            }
            else
            {
                threePrimePos = read.untrimmedStart();
                fivePrimePos = read.untrimmedEnd();
            }

            if(minFivePrimePos < 0 || fivePrimePos < minFivePrimePos)
                minFivePrimePos = fivePrimePos;

            if(minThreePrimePos < 0 || threePrimePos < minThreePrimePos)
                minThreePrimePos = threePrimePos;

            maxFivePrimePos = max(maxFivePrimePos, fivePrimePos);
            maxThreePrimePos = max(maxThreePrimePos, threePrimePos);
        }

        return (maxFivePrimePos - minFivePrimePos) + (maxThreePrimePos - minThreePrimePos) > SBX_PRIME_POSITION_RANGE_THRESHOLD;
    }
}
