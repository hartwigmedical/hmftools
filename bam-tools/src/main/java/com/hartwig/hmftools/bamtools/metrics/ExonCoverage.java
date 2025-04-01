package com.hartwig.hmftools.bamtools.metrics;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.VectorUtils;

public class ExonCoverage extends GeneRegion
{
    private final int[] mBaseCoverage;

    public ExonCoverage(final GeneRegion geneRegion)
    {
        super(geneRegion.Chromosome, geneRegion.start(), geneRegion.end(), geneRegion.GeneName, geneRegion.ExonRank);
        mBaseCoverage = new int[geneRegion.baseLength()];
    }

    public int[] coverage() { return mBaseCoverage;}

    public void processRead(int readStartPos, int readEndPos)
    {
        int startPosition = Math.max(start(), readStartPos);
        int endPosition = Math.min(end(), readEndPos);

        int startIndex = index(startPosition);
        int endIndex = index(endPosition);

        for(int i = startIndex; i <= endIndex; i++)
        {
            mBaseCoverage[i] += 1;
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

    public String toString() { return format("rank(%d) region(%s)", ExonRank, ((ChrBaseRegion)this)); }
}
