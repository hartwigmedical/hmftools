package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.WEAK_ASSEMBLY_UNPAIRED_LONG_EXT_FACTOR;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.readJunctionExtensionLength;
import static com.hartwig.hmftools.esvee.common.CommonUtils.isMediumBaseQual;
import static com.hartwig.hmftools.esvee.common.SvConstants.isUltima;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.read.ReadAdjustments;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;

public final class SeqTechUtils
{
    public static final int ULTIMA_READ_MISMATCH_MEDIUM_REPEAT_COUNT = 7;
    public static final int ULTIMA_READ_MISMATCH_LONG_REPEAT_COUNT = 11;

    public static final int SBX_PRIME_POSITION_RANGE_THRESHOLD = 6;
    public static final int SBX_MEDIUM_QUAL_DESYNC_COUNT = 150;

    public static final int UNPAIRED_READ_PROXIMATE_JUNCTION_DISTANCE = 100; // rather than BAM sampling for now

    public static void setSeqTechSpecifics()
    {
        if(isUltima())
        {
            AssemblyConstants.READ_MISMATCH_MEDIUM_REPEAT_COUNT = ULTIMA_READ_MISMATCH_MEDIUM_REPEAT_COUNT;
            AssemblyConstants.READ_MISMATCH_LONG_REPEAT_COUNT = ULTIMA_READ_MISMATCH_LONG_REPEAT_COUNT;
        }

        AssemblyConstants.PROXIMATE_JUNCTION_DISTANCE = UNPAIRED_READ_PROXIMATE_JUNCTION_DISTANCE;
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

    private static boolean withinLikelyDuplicateRange(final Read first, final Read second)
    {
        if(first.orientation() != second.orientation())
            return false;

        int firstStart = first.unclippedStart() - first.trimCountStart();
        int firstEnd = first.unclippedEnd() + first.trimCountEnd();
        int secondStart = second.unclippedStart() - second.trimCountStart();
        int secondEnd = second.unclippedEnd() + second.trimCountEnd();

        return abs(firstStart - secondStart) + abs(firstEnd - secondEnd) <= SBX_PRIME_POSITION_RANGE_THRESHOLD;
    }

    public static List<Read> findSbxPossibleDuplicates(final Junction junction, final List<Read> reads)
    {
        List<Integer> extensionLengths = Lists.newArrayListWithCapacity(reads.size());

        for(Read read : reads)
        {
            extensionLengths.add(readJunctionExtensionLength(read, junction));
        }

        // look for a group of reads with significantly longer extension lengths and likely duplicates of each other
        List<Integer> sortedLengths = Lists.newArrayList(extensionLengths);
        Collections.sort(sortedLengths, Collections.reverseOrder());

        int minLongLength = -1;
        int readsAboveLongLength = 0;
        for(int i = 0; i < sortedLengths.size() - 1; ++i)
        {
            int extLength = sortedLengths.get(i);
            int nextLength = sortedLengths.get(i + 1);

            if(extLength > nextLength * WEAK_ASSEMBLY_UNPAIRED_LONG_EXT_FACTOR)
            {
                minLongLength = extLength;
                readsAboveLongLength = i + 1;
                break;
            }
        }

        if(readsAboveLongLength < 2)
            return Collections.emptyList();

        List<Read> candidateDuplicates = Lists.newArrayListWithCapacity(readsAboveLongLength);

        for(int r = 0; r < reads.size(); ++r)
        {
            if(extensionLengths.get(r) >= minLongLength)
            {
                candidateDuplicates.add(reads.get(r));
            }
        }

        List<Read> likelyDuplicates = Lists.newArrayListWithCapacity(readsAboveLongLength - 1);

        for(int i = 0; i < candidateDuplicates.size() - 1; ++i)
        {
            Read read = candidateDuplicates.get(i);

            if(likelyDuplicates.contains(read))
                continue;

            for(int j = i + 1; j < candidateDuplicates.size(); ++j)
            {
                Read nextRead = candidateDuplicates.get(j);

                if(withinLikelyDuplicateRange(read, nextRead))
                {
                    likelyDuplicates.add(nextRead);
                }
            }
        }

        return likelyDuplicates;
    }

    public static int setSbxMediumQualCount(final Read read)
    {
        int mediumCount = 0;

        for(int i = 0; i < read.getBaseQuality().length; ++i)
        {
            if(isMediumBaseQual(read.getBaseQuality()[i]))
                ++mediumCount;
        }

        return mediumCount;
    }
}
