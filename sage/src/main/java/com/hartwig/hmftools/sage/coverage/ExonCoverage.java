package com.hartwig.hmftools.sage.coverage;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class ExonCoverage
{
    private final NamedBed mExon;
    private final int[] mBaseCoverage;

    ExonCoverage(final NamedBed exon)
    {
        mExon = exon;
        mBaseCoverage = new int[exon.bases()];
    }

    public int[] coverage()
    {
        return mBaseCoverage;
    }
    public String gene()
    {
        return mExon.name();
    }

    public String chromosome()
    {
        return mExon.chromosome();
    }

    public int start() { return mExon.start(); }

    public int end() { return mExon.end(); }

    public void processRead(int readStartPos, int readEndPos)
    {
        if(!positionsOverlap(readStartPos, readEndPos, mExon.start(), mExon.end()))
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


    private int index(int position) { return position - start(); }
}
