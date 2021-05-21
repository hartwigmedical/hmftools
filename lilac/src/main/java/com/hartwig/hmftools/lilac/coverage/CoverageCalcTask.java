package com.hartwig.hmftools.lilac.coverage;

import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

public class CoverageCalcTask implements Callable
{
    private final List<FragmentAlleles> mFragmentAlleles;
    private final List<HlaComplex> mComplexes;
    private List<HlaComplexCoverage> mCoverageResults;

    private final Map<HlaAllele,List<FragmentAlleles>> mAlleleFragmentMap;

    private final PerformanceCounter mPerfCounterFilter;
    private final PerformanceCounter mPerfCounterCoverage;

    public CoverageCalcTask(
            final List<FragmentAlleles> fragmentAlleles, final List<HlaComplex> complexes,
            Map<HlaAllele,List<FragmentAlleles>> alleleFragmentMap)
    {
        mFragmentAlleles = fragmentAlleles;
        mComplexes = complexes;
        mCoverageResults = Lists.newArrayList();

        mAlleleFragmentMap = alleleFragmentMap;

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
        mComplexes.forEach(x -> mCoverageResults.add(proteinCoverage(mFragmentAlleles, x.getAlleles())));
        return (long)0;
    }

    private HlaComplexCoverage proteinCoverage(final List<FragmentAlleles> fragmentAlleles, final List<HlaAllele> alleles)
    {
        // gather all fragments which cover this set of alleles
        mPerfCounterFilter.start();
        List<FragmentAlleles> filteredFragments = FragmentAlleles.filter(fragmentAlleles, alleles);
        mPerfCounterFilter.stop();

        // tally up protein-supporting coverage into the complex coverage counts
        mPerfCounterCoverage.start();
        HlaComplexCoverage coverage = HlaComplexCoverage.create(HlaAlleleCoverage.proteinCoverage(filteredFragments));
        mPerfCounterCoverage.stop();

        return coverage;
    }

    public void logPerfResults()
    {
        /*
        LL_LOGGER.debug(String.format("filter perf: count(%d) avg(%.4f) max(%.4f)",
                mPerfCounterFilter.getSampleCount(), mPerfCounterFilter.getAvgTime(), mPerfCounterFilter.getMaxTime()));

        LL_LOGGER.debug(String.format("coverage perf: count(%d) avg(%.4f) max(%.4f)",
                mPerfCounterCoverage.getSampleCount(), mPerfCounterCoverage.getAvgTime(), mPerfCounterCoverage.getMaxTime()));
        */
    }

}
