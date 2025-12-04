package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.common.SvConstants.isUltima;

import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.read.ReadAdjustments;

public final class SeqTechUtils
{
    public static final int READ_MISMATCH_MEDIUM_REPEAT_COUNT_ULTIMA = 7;
    public static final int READ_MISMATCH_LONG_REPEAT_COUNT_ULTIMA = 11;

    public static void setSeqTechSpecifics()
    {
        if(isUltima())
        {
            AssemblyConstants.READ_MISMATCH_MEDIUM_REPEAT_COUNT = READ_MISMATCH_MEDIUM_REPEAT_COUNT_ULTIMA;
            AssemblyConstants.READ_MISMATCH_LONG_REPEAT_COUNT = READ_MISMATCH_LONG_REPEAT_COUNT_ULTIMA;
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
}
