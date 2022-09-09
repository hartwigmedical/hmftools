package com.hartwig.hmftools.lilac.coverage;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

public class CoverageCalcTask implements Callable<Long>
{
    private final int mId;
    private final List<HlaComplex> mComplexes;
    private List<ComplexCoverage> mCoverageResults;

    private final FragmentAlleleMatrix mFragAlleleMatrix;
    private final double mTopScorePercDiff;
    private int mMaxFragments;
    private int mLowScoreCount;

    private final PerformanceCounter mPerfCounter;

    // private static final int CULL_COMPLEX_COUNT = 500000;
    private static final int CULL_COMPLEX_COUNT = 100000;
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
        mPerfCounter = new PerformanceCounter("CalcCoverage");
    }

    public ComplexCoverage getCoverage() { return mCoverageResults.get(0); }
    public List<ComplexCoverage> getCoverageResults() { return mCoverageResults; }
    public PerformanceCounter getPerfCounter() { return mPerfCounter; }

    @Override
    public Long call()
    {
        boolean checkCull = mComplexes.size() >= CULL_COMPLEX_COUNT;

        mPerfCounter.start();

        for(int i = 0; i < mComplexes.size(); ++i)
        {
            if(checkCull && i > 0 && (i % CULL_COMPLEX_COUNT) == 0)
            {
                LL_LOGGER.debug(String.format("thread %d: complexes(%d) processed, discard(%d, %.2f%%)",
                        mId, i, mLowScoreCount, 100.0 * mLowScoreCount / i));

                mPerfCounter.stop();
                System.gc();
                mPerfCounter.start();
            }

            List<AlleleCoverage> alleleCoverage = mFragAlleleMatrix.create(mComplexes.get(i));
            int totalFragments = calcTotalFragments(alleleCoverage);

            if(checkCull && canCull(totalFragments)) // result.TotalCoverage
                continue;

            ComplexCoverage result = ComplexCoverage.create(alleleCoverage);

            mCoverageResults.add(result);
        }

        mPerfCounter.stop();

        return (long)0;
    }

    private int calcTotalFragments(List<AlleleCoverage> alleleCoverage)
    {
        int unique = 0;
        double shared = 0.0;
        double wild = 0.0;

        for(AlleleCoverage coverage : alleleCoverage)
        {
            unique += coverage.UniqueCoverage;
            shared += coverage.SharedCoverage;
            wild += coverage.WildCoverage;
        }

        return unique + (int)round(shared) + (int)round(wild);
    }

    private boolean canCull(final int totalCoverage)
    {
        if(totalCoverage > mMaxFragments)
        {
            mMaxFragments = totalCoverage;
            return false;
        }

        if(mMaxFragments - totalCoverage < MIN_FRAG_DIFF)
            return false;

        if(totalCoverage > mMaxFragments * (1 - mTopScorePercDiff))
            return false;

        ++mLowScoreCount;
        return true;
    }

    private ComplexCoverage calcCoverage(final HlaComplex complex)
    {
        List<AlleleCoverage> alleleCoverage = mFragAlleleMatrix.create(complex);
        return ComplexCoverage.create(alleleCoverage);
    }
}
