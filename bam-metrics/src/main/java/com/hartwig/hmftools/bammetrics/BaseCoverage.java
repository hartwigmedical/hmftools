package com.hartwig.hmftools.bammetrics;

import static com.hartwig.hmftools.bammetrics.FilterType.MAX_COVERAGE;
import static com.hartwig.hmftools.bammetrics.FilterType.OVERLAPPED;
import static com.hartwig.hmftools.bammetrics.FilterType.UNFILTERED;
import static com.hartwig.hmftools.bammetrics.FilterType.LOW_BASE_QUAL;
import static com.hartwig.hmftools.bammetrics.FilterType.LOW_MAP_QUAL;
import static com.hartwig.hmftools.bammetrics.FilterType.DUPLICATE;
import static com.hartwig.hmftools.bammetrics.FilterType.MATE_UNMAPPED;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.mateUnmapped;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import static htsjdk.samtools.CigarOperator.M;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class BaseCoverage
{
    private final BmConfig mConfig;
    private final int mRegionSize;
    private int mRegionStart;

    private final int[] mBaseDepth;
    private final long[] mFilterTypeCounts;

    public BaseCoverage(final BmConfig config, int regionStart, int regionEnd)
    {
        mConfig = config;
        mRegionSize = regionEnd - regionStart + 1;
        mRegionStart = regionStart;
        mBaseDepth = new int[mRegionSize];
        mFilterTypeCounts = new long[FilterType.values().length];
    }

    public void processRead(final SAMRecord record, final int[] mateBaseCoords)
    {
        // some filters exclude all matched bases
        int alignedBases = record.getCigar().getCigarElements().stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();

        if(record.getMappingQuality() < mConfig.MapQualityThreshold)
        {
            mFilterTypeCounts[LOW_MAP_QUAL.ordinal()] += alignedBases;
            return;
        }

        if(record.getDuplicateReadFlag())
        {
            mFilterTypeCounts[DUPLICATE.ordinal()] += alignedBases;
            return;
        }

        if(mateUnmapped(record))
        {
            mFilterTypeCounts[MATE_UNMAPPED.ordinal()] += alignedBases;
            return;
        }

        int position = record.getAlignmentStart();
        int readIndex = 0;

        for(CigarElement element : record.getCigar().getCigarElements())
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
                            record, position, readIndex, element.getLength(), mateBaseCoords);
                    break;

                default:
                    break;
            }
        }
    }

    private void processMatchedBases(
            final SAMRecord record, int posStart, int readIndexStart, int matchLength, final int[] mateBaseCoords)
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

            boolean lowBaseQual = record.getBaseQualities()[readIndex] < mConfig.BaseQualityThreshold;

            boolean overlapped = mateBaseCoords != null && positionWithin(position, mateBaseCoords[SE_START], mateBaseCoords[SE_END]);

            boolean exceedsCoverage = mBaseDepth[baseIndex] >= mConfig.MaxCoverage;

            boolean passFilters = !lowBaseQual && !exceedsCoverage && !overlapped;

            if(passFilters)
            {
                ++mBaseDepth[baseIndex];
                ++mFilterTypeCounts[UNFILTERED.ordinal()];
            }
            else
            {
                if(lowBaseQual)
                    ++mFilterTypeCounts[LOW_BASE_QUAL.ordinal()];
                else if(overlapped)
                    ++mFilterTypeCounts[OVERLAPPED.ordinal()];
                else if(exceedsCoverage)
                    ++mFilterTypeCounts[MAX_COVERAGE.ordinal()];
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
}
