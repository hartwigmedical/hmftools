package com.hartwig.hmftools.lilac.coverage;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;

public class CoverageCalcTask implements Callable
{
    private final int mId;
    private final List<HlaComplex> mComplexes;
    private List<ComplexCoverage> mCoverageResults;

    private final FragmentAlleleMatrix mFragAlleleMatrix;
    private final double mTopScorePercDiff;
    private int mMaxFragments;
    private int mLowScoreCount;

    private static final int CULL_COMPLEX_COUNT = 500000;
    private static final int LOG_COMPLEX_COUNT = 250000;
    private static final int MIN_FRAG_DIFF = 40;

    public CoverageCalcTask(final int id, final List<HlaComplex> complexes, FragmentAlleleMatrix fragAlleleMatrix, double topScoreThreshold)
    {
        mId = id;
        mComplexes = complexes;
        mCoverageResults = Lists.newArrayList();

        mFragAlleleMatrix = fragAlleleMatrix;
        mTopScorePercDiff = min(topScoreThreshold * 5, 0.99);
        mMaxFragments = 0;
        mLowScoreCount = 0;
    }

    public ComplexCoverage getCoverage() { return mCoverageResults.get(0); }
    public List<ComplexCoverage> getCoverageResults() { return mCoverageResults; }

    @Override
    public Long call()
    {
        boolean checkCull = mComplexes.size() >= CULL_COMPLEX_COUNT;

        for(int i = 0; i < mComplexes.size(); ++i)
        {
            if(checkCull && i > 0 && (i % LOG_COMPLEX_COUNT) == 0)
            {
                LL_LOGGER.debug(String.format("thread %d: complexes(%d) processed, discard(%d, %.0f%%)",
                        mId, i, mLowScoreCount, 100.0 * mLowScoreCount / i));
            }

            ComplexCoverage result = calcCoverage(mComplexes.get(i));

            if(checkCull && canCull(result))
                continue;

            mCoverageResults.add(result);
        }

        return (long)0;
    }

    private boolean canCull(final ComplexCoverage result)
    {
        if(result.TotalCoverage > mMaxFragments)
        {
            mMaxFragments = result.TotalCoverage;
            return false;
        }

        if(mMaxFragments - result.TotalCoverage < MIN_FRAG_DIFF)
            return false;

        if(result.TotalCoverage > mMaxFragments * (1 - mTopScorePercDiff))
            return false;

        ++mLowScoreCount;
        return true;
    }

    private ComplexCoverage calcCoverage(final HlaComplex complex)
    {
        List<AlleleCoverage> alleleCoverage = mFragAlleleMatrix.create(complex);
        return ComplexCoverage.create(alleleCoverage);
    }

    public void logPerfResults()
    {
        /*
        LL_LOGGER.debug(String.format("filter perf: count(%d) avg(%.4f) max(%.4f)",
                mPerfCounterFilter.getSampleCount(), mPerfCounterFilter.getAvgTime(), mPerfCounterFilter.getMaxTime()));

        LL_LOGGER.debug(String.format("coverage perf: count(%d) avg(%.6f) max(%.6f)",
                mPerfCounterCoverage.getSampleCount(), mPerfCounterCoverage.getAvgTime(), mPerfCounterCoverage.getMaxTime()));
        */
    }

}
