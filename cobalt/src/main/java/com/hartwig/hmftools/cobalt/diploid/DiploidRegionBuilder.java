package com.hartwig.hmftools.cobalt.diploid;

import static com.hartwig.hmftools.cobalt.CobaltConstants.WINDOW_SIZE;

import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public class DiploidRegionBuilder implements Consumer<DiploidCount>
{
    private final List<GenomeRegion> mResult;

    private final double mCutoff;
    private final int mMaleSamples;
    private final int mFemaleSamples;
    private long mTotalDiploidBases;

    private String mChromosome = "";
    private long mStart = 0;
    private long mEnd = 0;
    private boolean mIsDiploid;

    public DiploidRegionBuilder(final double cutoff, final int femaleSamples, final int maleSamples)
    {
        mResult = Lists.newArrayList();

        mCutoff = cutoff;
        mMaleSamples = maleSamples;
        mFemaleSamples = femaleSamples;
        mTotalDiploidBases = 0;
        mStart = 0;
        mEnd = 0;
        mChromosome = "";
    }

    @Override
    public void accept(@NotNull DiploidCount count)
    {
        int samples = count.chromosome().equals("Y") || count.chromosome().equals("chrY") ? mMaleSamples : mFemaleSamples;
        boolean isCountDiploid = Doubles.greaterOrEqual(count.proportionIsDiploid(samples), mCutoff);
        if(!count.chromosome().equals(mChromosome) || isCountDiploid != mIsDiploid)
        {
            finaliseCurrent();
            mChromosome = count.chromosome();
            mStart = count.position();
            mIsDiploid = isCountDiploid;
        }

        mEnd = count.position() + WINDOW_SIZE - 1;
    }

    private void finaliseCurrent()
    {
        if(mIsDiploid)
        {
            mResult.add(GenomeRegions.create(mChromosome, mStart, mEnd));
            mTotalDiploidBases += mEnd - mStart + 1;
        }
        mIsDiploid = false;
    }

    public long getTotalDiploidBases()
    {
        return mTotalDiploidBases;
    }

    public List<GenomeRegion> build()
    {
        finaliseCurrent();
        return mResult;
    }
}
