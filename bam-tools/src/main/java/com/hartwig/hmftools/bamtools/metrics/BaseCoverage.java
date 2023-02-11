package com.hartwig.hmftools.bamtools.metrics;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import static htsjdk.samtools.CigarOperator.M;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class BaseCoverage
{
    private final MetricsConfig mConfig;
    private final int mRegionSize;
    private int mRegionStart;

    private final int[] mBaseDepth;
    private final long[] mFilterTypeCounts;

    public BaseCoverage(final MetricsConfig config, int regionStart, int regionEnd)
    {
        mConfig = config;
        mRegionSize = regionEnd - regionStart + 1;
        mRegionStart = regionStart;
        mBaseDepth = new int[mRegionSize];
        mFilterTypeCounts = new long[FilterType.values().length];
    }

    public void processRead(final SAMRecord read, final List<int[]> mateBaseCoords)
    {
        // some filters exclude all matched bases
        int alignedBases = read.getCigar().getCigarElements().stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();

        if(read.getMappingQuality() < mConfig.MapQualityThreshold)
        {
            mFilterTypeCounts[FilterType.LOW_MAP_QUAL.ordinal()] += alignedBases;
            return;
        }

        if(read.getDuplicateReadFlag())
        {
            mFilterTypeCounts[FilterType.DUPLICATE.ordinal()] += alignedBases;
            return;
        }

        if(!read.getReadPairedFlag() || read.getMateUnmappedFlag())
        {
            mFilterTypeCounts[FilterType.MATE_UNMAPPED.ordinal()] += alignedBases;
            return;
        }

        int position = read.getAlignmentStart();
        int readIndex = 0;

        for(CigarElement element : read.getCigar().getCigarElements())
        {
            switch(element.getOperator())
            {
                case S:
                case H:
                    readIndex += element.getLength();
                    break;

                case D:
                case N:
                    position += element.getLength();
                    break;

                case I:
                    readIndex += element.getLength();
                    break;

                case M:
                    processMatchedBases(
                            read, position, readIndex, element.getLength(), mateBaseCoords);
                    position += element.getLength();
                    readIndex += element.getLength();
                    break;

                default:
                    break;
            }
        }
    }

    private void processMatchedBases(
            final SAMRecord read, int posStart, int readIndexStart, int matchLength, final List<int[]> mateBaseCoords)
    {
        for(int i = 0; i < matchLength; ++i)
        {
            int position = posStart + i;

            if(position < mRegionStart)
                continue;

            if(position >= mRegionStart + mRegionSize)
                break;

            int readIndex = readIndexStart + i;
            int baseIndex = position - mRegionStart;

            boolean lowBaseQual = read.getBaseQualities()[readIndex] < mConfig.BaseQualityThreshold;

            boolean overlapped = mateBaseCoords != null && mateBaseCoords.stream().anyMatch(x -> positionWithin(position, x[SE_START], x[SE_END]));

            boolean exceedsCoverage = mBaseDepth[baseIndex] >= mConfig.MaxCoverage;

            boolean passFilters = !lowBaseQual && !exceedsCoverage && !overlapped;

            if(passFilters)
            {
                ++mBaseDepth[baseIndex];
                ++mFilterTypeCounts[FilterType.UNFILTERED.ordinal()];
            }
            else
            {
                if(lowBaseQual)
                    ++mFilterTypeCounts[FilterType.LOW_BASE_QUAL.ordinal()];
                else if(overlapped)
                    ++mFilterTypeCounts[FilterType.OVERLAPPED.ordinal()];
                else if(exceedsCoverage)
                    ++mFilterTypeCounts[FilterType.MAX_COVERAGE.ordinal()];
            }
        }
    }

    public Metrics createMetrics()
    {
        Metrics metrics = new Metrics(mConfig.MaxCoverage);

        long coverageBases = 0;

        for(int i = 0; i < mBaseDepth.length; ++i)
        {
            int coverage = mBaseDepth[i];

            if(coverage == 0)
            {
                if(mConfig.ExcludeZeroCoverage)
                    continue;
            }
            else
            {
                ++coverageBases;
            }

            ++metrics.CoverageFrequency[coverage];
        }

        for(FilterType type : FilterType.values())
        {
            metrics.FilterTypeCounts[type.ordinal()] += mFilterTypeCounts[type.ordinal()];
        }

        metrics.addCoverageBases(coverageBases);

        return metrics;
    }

    public void clear()
    {
        for(int i = 0; i < mBaseDepth.length; ++i)
        {
            mBaseDepth[i] = 0;
        }

        for(int i = 0; i < mFilterTypeCounts.length; ++i)
        {
            mFilterTypeCounts[i] = 0;
        }
    }

    @VisibleForTesting
    public int[] baseDepth() { return mBaseDepth; }
}
