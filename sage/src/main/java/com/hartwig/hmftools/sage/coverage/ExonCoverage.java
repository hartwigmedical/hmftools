package com.hartwig.hmftools.sage.coverage;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

import org.jetbrains.annotations.NotNull;

public class ExonCoverage implements Consumer<BaseRegion>
{
    private final NamedBed mExon;
    private final int[] mBaseCoverage;

    ExonCoverage(final NamedBed exon)
    {
        mExon = exon;
        mBaseCoverage = new int[(int) exon.bases()];
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

    public long start()
    {
        return mExon.start();
    }

    public long end()
    {
        return mExon.end();
    }

    @Override
    public void accept(final BaseRegion alignment)
    {
        if(alignment.start() <= mExon.end() && alignment.end() >= mExon.start())
        {
            int startPosition = (int) Math.max(start(), alignment.start());
            int endPosition = (int) Math.min(end(), alignment.end());

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
    }


    private int index(int position)
    {
        return (int) (position - start());
    }
}
