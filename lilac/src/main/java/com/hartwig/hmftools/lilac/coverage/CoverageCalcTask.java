package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.common.sigs.DataUtils.doublesEqual;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

public class CoverageCalcTask implements Callable
{
    private final List<FragmentAlleles> mFragmentAlleles;
    private final List<HlaComplex> mComplexes;
    private List<HlaComplexCoverage> mCoverageResults;

    private final FragmentAlleleMatrix mFragAlleleMatrix;

    private final PerformanceCounter mPerfCounterFilter;
    private final PerformanceCounter mPerfCounterCoverage;

    public CoverageCalcTask(
            final List<FragmentAlleles> fragmentAlleles, final List<HlaComplex> complexes,
            FragmentAlleleMatrix fragAlleleMatrix)
    {
        mFragmentAlleles = fragmentAlleles;
        mComplexes = complexes;
        mCoverageResults = Lists.newArrayList();

        mFragAlleleMatrix = fragAlleleMatrix;

        mPerfCounterFilter = new PerformanceCounter("ComplexFilter");
        mPerfCounterCoverage = new PerformanceCounter("ComplexCoverage");
    }

    public CoverageCalcTask(final List<FragmentAlleles> fragmentAlleles, final HlaComplex complex)
    {
        this(fragmentAlleles, Lists.newArrayList(complex), null);
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

        boolean compareWithOld = false;

        if(compareWithOld)
        {
            // mPerfCounterCoverage.start();
            // mPerfCounterCoverage.stop();

            // gather all fragments which cover this set of alleles
            List<FragmentAlleles> filteredFragments = FragmentAlleles.filter(mFragmentAlleles, complex.getAlleles());
            List<HlaAlleleCoverage> origCoverage = HlaAlleleCoverage.proteinCoverage(filteredFragments);

            for(int i = 0; i < origCoverage.size(); ++i)
            {
                HlaAlleleCoverage origCov = origCoverage.get(i);
                HlaAlleleCoverage newCov = alleleCoverage.stream().filter(x -> x.Allele.equals(origCov.Allele)).findFirst().orElse(null);

                if(origCov.UniqueCoverage != newCov.UniqueCoverage
                || !doublesEqual(origCov.TotalCoverage, newCov.TotalCoverage))
                {
                    LL_LOGGER.warn("differing cov");
                }
            }

            return HlaComplexCoverage.create(origCoverage);
        }

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
