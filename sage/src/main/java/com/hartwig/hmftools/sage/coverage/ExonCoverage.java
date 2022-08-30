package com.hartwig.hmftools.sage.coverage;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.utils.VectorUtils;

import org.apache.commons.compress.utils.Lists;

public class ExonCoverage
{
    private final GenomeRegion mRegion;
    private final int mExonRank;
    private final int[] mBaseCoverage;

    public ExonCoverage(final GenomeRegion region, final int exonRank)
    {
        mRegion = region;
        mExonRank = exonRank;
        mBaseCoverage = new int[region.bases()];
    }

    public final GenomeRegion region() { return mRegion; }
    public int exonRank() { return mExonRank; }

    public int[] coverage() { return mBaseCoverage;}

    public String chromosome() { return mRegion.chromosome(); }
    public int start() { return mRegion.start(); }
    public int end() { return mRegion.end(); }

    public void processRead(int readStartPos, int readEndPos)
    {
        if(!positionsOverlap(readStartPos, readEndPos, mRegion.start(), mRegion.end()))
            return;

        int startPosition = Math.max(start(), readStartPos);
        int endPosition = Math.min(end(), readEndPos);

        int startIndex = index(startPosition);
        int endIndex = index(endPosition);

        synchronized (mBaseCoverage)
        {
            for(int i = startIndex; i <= endIndex; i++)
            {
                mBaseCoverage[i] += 1;
            }
        }
    }

    public double medianDepth()
    {
        List<Double> depths = Lists.newArrayList();

        for(int i = 0; i < mBaseCoverage.length; ++i)
        {
            VectorUtils.optimisedAdd(depths, mBaseCoverage[i], true);
        }

        int medianIndex = mBaseCoverage.length / 2;
        return depths.get(medianIndex);
    }

    private int index(int position) { return position - start(); }

    public String toString() { return format("rank(%d) region(%s)", mExonRank, mRegion.toString()); }
}
