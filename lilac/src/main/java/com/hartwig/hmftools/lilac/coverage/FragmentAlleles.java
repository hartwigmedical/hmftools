package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_A;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_Y_FRAGMENT_THRESHOLD;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.CANDIDATE;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.HLA_Y;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.NO_HET_LOCI;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.UNMATCHED_AMINO_ACID;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceLoci.filterExonBoundaryWildcards;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceLoci.filterWildcards;
import static com.hartwig.hmftools.lilac.seq.SequenceMatchType.FULL;
import static com.hartwig.hmftools.lilac.seq.SequenceMatchType.WILD;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.SequenceCount;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.fragment.FragmentScope;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.seq.SequenceMatchType;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
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
