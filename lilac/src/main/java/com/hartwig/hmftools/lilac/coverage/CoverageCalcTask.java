package com.hartwig.hmftools.lilac.coverage;

import java.util.List;
import java.util.concurrent.Callable;

import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.read.FragmentAlleles;

public class CoverageCalcTask implements Callable
{
    private final List<FragmentAlleles> mFragmentAlleles;
    private final HlaComplex mComplex;
    private HlaComplexCoverage mCoverage;

    public CoverageCalcTask(final List<FragmentAlleles> fragmentAlleles, final HlaComplex complex)
    {
        mFragmentAlleles = fragmentAlleles;
        mComplex = complex;
        mCoverage = null;
    }

    public HlaComplexCoverage getCoverage() { return mCoverage; }


    @Override
    public Long call()
    {
        mCoverage = proteinCoverage(mFragmentAlleles, mComplex.getAlleles());
        return (long)0;
    }

    public static HlaComplexCoverage proteinCoverage(final List<FragmentAlleles> fragmentAlleles, final List<HlaAllele> alleles)
    {
        List<FragmentAlleles> filteredFragments = fragmentAlleles(fragmentAlleles, alleles);
        return HlaComplexCoverage.create(HlaAlleleCoverage.proteinCoverage(filteredFragments));
    }

    public static List<FragmentAlleles> fragmentAlleles(final List<FragmentAlleles> fragmentAlleles, final List<HlaAllele> alleles)
    {
        return FragmentAlleles.filter(fragmentAlleles, alleles);
    }

    public static HlaComplexCoverage groupCoverage(final List<FragmentAlleles> fragmentAlleles, final List<HlaAllele> alleles)
    {
        List<FragmentAlleles> filteredFragments = fragmentAlleles(fragmentAlleles, alleles);
        return HlaComplexCoverage.create(HlaAlleleCoverage.groupCoverage(filteredFragments));
    }


}
