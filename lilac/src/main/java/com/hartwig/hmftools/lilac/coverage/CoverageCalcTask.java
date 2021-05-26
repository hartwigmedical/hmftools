package com.hartwig.hmftools.lilac.coverage;

import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

public class CoverageCalcTask implements Callable
{
    private final List<HlaComplex> mComplexes;
    private List<HlaComplexCoverage> mCoverageResults;

    private final FragmentAlleleMatrix mFragAlleleMatrix;

    private final PerformanceCounter mPerfCounter;

    public CoverageCalcTask(final List<HlaComplex> complexes, FragmentAlleleMatrix fragAlleleMatrix)
    {
        mComplexes = complexes;
        mCoverageResults = Lists.newArrayList();

        mFragAlleleMatrix = fragAlleleMatrix;

        mPerfCounter = new PerformanceCounter("ComplexCalcs");
    }

    public HlaComplexCoverage getCoverage() { return mCoverageResults.get(0); }
    public List<HlaComplexCoverage> getCoverageResults() { return mCoverageResults; }

    @Override
    public Long call()
    {
        mComplexes.forEach(x -> mCoverageResults.add(calcCoverage(x)));
        return (long)0;
    }

    private HlaComplexCoverage calcCoverage(final HlaComplex complex)
    {
        List<HlaAlleleCoverage> alleleCoverage = mFragAlleleMatrix.create(complex);
        return HlaComplexCoverage.create(alleleCoverage);
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
