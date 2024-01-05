package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.common.AssemblySequence;
import com.hartwig.hmftools.esvee.common.Direction;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.old.PrimaryAssembly;
import com.hartwig.hmftools.esvee.read.ReadFilters;

public class PrimaryAssembler
{
    private final SvConfig mConfig;

    private final Junction mJunction;

    private int mNextAssemblyNumber = 1;

    public PrimaryAssembler(final SvConfig config, final Junction junction)
    {
        mConfig = config;
        mJunction = junction;
    }

    public List<PrimaryAssembly> processJunction(final List<Read> rawReads)
    {
        // FIXME:
        final List<Read> realignedReads = rawReads.stream()
                // .map(alignment -> realignForJunction(alignment, mJunction))
                .collect(Collectors.toList());

        final List<Read> withLowQAlignments = realignedReads.stream()
                .filter(alignment -> ReadFilters.recordSoftClipsNearJunction(alignment, mJunction)) // mCounters.ReadsSoftClippedAtJunction
                .collect(Collectors.toList());

        final List<Read> filteredAlignments = withLowQAlignments.stream()
                .filter(alignment -> ReadFilters.isRecordAverageQualityPastJunctionAbove(alignment, mJunction, SvConstants.AVG_BASE_QUAL_THRESHOLD)) // mCounters.ReadsPassingJunctionQualityThreshold
                .filter(alignment -> ReadFilters.hasAcceptableMapQ(alignment, SvConstants.READ_FILTER_MIN_JUNCTION_MAPQ)) // mCounters.HasAcceptableMapQ
                .filter(ReadFilters::isNotBadlyMapped) // mCounters.WellMapped
                .collect(Collectors.toList());

        if(filteredAlignments.isEmpty())
            return List.of(); // There are no reads of acceptable quality supporting this junction

        final List<PrimaryAssembly> initialAssemblies = createInitialAssemblies(filteredAlignments);

        return Collections.emptyList();
    }

    private String nextAssemblyName()
    {
        return String.format("%s:%s%s:%s", mJunction.Chromosome, mJunction.Position,
                mJunction.direction() == Direction.FORWARDS ? "F" : "R", mNextAssemblyNumber++);
    }

    private List<PrimaryAssembly> createInitialAssemblies(final List<Read> junctionReads)
    {
        int minAlignedPosition = mJunction.Position;
        int maxAlignedPosition = mJunction.Position;

        Read maxJunctionBaseQualRead = null;
        int maxJunctionBaseQualTotal = 0;

        for(Read read : junctionReads)
        {
            if(mJunction.direction() == Direction.FORWARDS)
            {
                maxAlignedPosition = max(maxAlignedPosition, read.getUnclippedEnd());
            }
            else
            {
                minAlignedPosition = min(minAlignedPosition, read.getUnclippedStart());
            }

            int junctionBaseQualTotal = readQualFromJunction(read, mJunction);

            if(junctionBaseQualTotal > maxJunctionBaseQualTotal)
            {
                maxJunctionBaseQualTotal = junctionBaseQualTotal;
                maxJunctionBaseQualRead = read;
            }
        }

        AssemblySequence assemblySequence = new AssemblySequence(mJunction, maxJunctionBaseQualRead, minAlignedPosition,  maxAlignedPosition);

        for(Read read : junctionReads)
        {
            if(read == maxJunctionBaseQualRead)
                continue;

            assemblySequence.tryAddRead(read);
        }

        return Collections.emptyList();
    }

    private static int readQualFromJunction(final Read read, final Junction junction)
    {
        int junctionReadIndex = read.getReadIndexAtReferencePosition(junction.Position, true);

        int readIndexStart;
        int readIndexEnd;

        if(junction.direction() == Direction.FORWARDS)
        {
            readIndexStart = junctionReadIndex;
            readIndexEnd = read.getLength() - 1;
        }
        else
        {
            readIndexStart = 0;
            readIndexEnd = junctionReadIndex;
        }

        int baseQualTotal = 0;

        for(int i = readIndexStart; i <= readIndexEnd; ++i)
        {
            baseQualTotal += read.getBaseQuality()[i];
        }

        return baseQualTotal;
    }
}
