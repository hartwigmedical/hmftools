package com.hartwig.hmftools.lilac.coverage;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class FragmentAlleles
{
    private final Fragment mFragment;
    private final Set<HlaAllele> mFull;
    private final Set<HlaAllele> mWild;

    public FragmentAlleles(final Fragment fragment, final List<HlaAllele> full, final List<HlaAllele> wild)
    {
        mFragment = fragment;
        mFull = Sets.newHashSet(full);
        mWild = Sets.newHashSet(wild);
    }

    public boolean contains(final HlaAllele allele)
    {
        return mFull.contains(allele) || mWild.contains(allele);
    }

    public final Fragment getFragment() { return mFragment; }

    public final Set<HlaAllele> getFull() { return mFull; }
    public final Set<HlaAllele> getWild() { return mWild; }

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
