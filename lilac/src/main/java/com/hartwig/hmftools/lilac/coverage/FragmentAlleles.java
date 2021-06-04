package com.hartwig.hmftools.lilac.coverage;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import java.util.List;
import java.util.stream.Collectors;

public class FragmentAlleles
{
    private final Fragment mFragment;
    private final List<HlaAllele> mFull;
    private final List<HlaAllele> mWild;

    public FragmentAlleles(
            final Fragment fragment, final List<HlaAllele> full, final List<HlaAllele> wild)
    {
        mFragment = fragment;
        mFull = full;
        mWild = wild;
    }

    public boolean contains(final HlaAllele allele)
    {
        return mFull.contains(allele) || mWild.contains(allele);
    }

    public final Fragment getFragment() { return mFragment; }

    public final List<HlaAllele> getFull() { return mFull; }
    public final List<HlaAllele> getWild() { return mWild; }

    public static List<FragmentAlleles> filter(final List<FragmentAlleles> fragAlleleList, final List<HlaAllele> alleles)
    {
        // gather any fragment allele which contains at least one of the specified alleles in its full or wild list,
        // then collecting any matching alleles in each of the three groups
        List<FragmentAlleles> matchedFragAlleles = Lists.newArrayList();

        for(FragmentAlleles fragAllele : fragAlleleList)
        {
            if(alleles.stream().anyMatch(x -> fragAllele.contains(x)))
            {
                matchedFragAlleles.add(new FragmentAlleles(
                        fragAllele.getFragment(),
                        alleles.stream().filter(x -> fragAllele.getFull().contains(x)).collect(Collectors.toList()),
                        alleles.stream().filter(x -> fragAllele.getWild().contains(x)).collect(Collectors.toList())));
            }
        }

        return matchedFragAlleles;
    }
}
