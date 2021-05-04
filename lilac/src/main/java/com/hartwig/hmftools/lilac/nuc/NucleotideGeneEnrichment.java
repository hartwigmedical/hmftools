package com.hartwig.hmftools.lilac.nuc;

import static com.hartwig.hmftools.lilac.LilacConstants.HLA_A;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_B;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_C;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;

public class NucleotideGeneEnrichment
{
    private final int mAbMinBoundary;
    private final int mAcMinBoundary;
    private final int mBcMinBoundary;

    public NucleotideGeneEnrichment(final Set<Integer> aBoundaries, final Set<Integer> bBoundaries, final Set<Integer> cBoundaries)
    {
        Set<Integer> abUnique = Sets.newHashSet();
        Set<Integer> acUnique = Sets.newHashSet();
        Set<Integer> bcUnique = Sets.newHashSet();

        aBoundaries.stream().filter(x -> !bBoundaries.contains(x)).forEach(x -> abUnique.add(x));
        bBoundaries.stream().filter(x -> !aBoundaries.contains(x)).forEach(x -> abUnique.add(x));
        mAbMinBoundary = abUnique.stream().mapToInt(x -> x).min().orElse(0);

        aBoundaries.stream().filter(x -> !bBoundaries.contains(x)).forEach(x -> acUnique.add(x));
        cBoundaries.stream().filter(x -> !aBoundaries.contains(x)).forEach(x -> acUnique.add(x));
        mAcMinBoundary = acUnique.stream().mapToInt(x -> x).min().orElse(0);

        bBoundaries.stream().filter(x -> !cBoundaries.contains(x)).forEach(x -> bcUnique.add(x));
        cBoundaries.stream().filter(x -> !bBoundaries.contains(x)).forEach(x -> bcUnique.add(x));
        mBcMinBoundary = bcUnique.stream().mapToInt(x -> x).min().orElse(0);
    }

    public final int getAFilterB() { return mAbMinBoundary; }
    public final int getAFilterC()
    {
        return mAcMinBoundary;
    }
    public final int getBFilterA()
    {
        return mAbMinBoundary;
    }
    public final int getBFilterC()
    {
        return mBcMinBoundary;
    }
    public final int getCFilterA()
    {
        return mAcMinBoundary;
    }
    public final int getCFilterB()
    {
        return mBcMinBoundary;
    }

    public final List<NucleotideFragment> enrich(final List<NucleotideFragment> fragments)
    {
        return fragments.stream().map(x -> enrich(x)).collect(Collectors.toList());
    }

    public final NucleotideFragment enrichGenes(final NucleotideFragment fragment)
    {
        return enrich(fragment);
    }

    public final NucleotideFragment enrich(final NucleotideFragment fragment)
    {
        if(fragment.containsIndel())
            return fragment;

        Set<String> genes = Sets.newHashSet();

        if(matchToA(fragment))
            genes.add(HLA_A);

        if(matchToB(fragment))
            genes.add(HLA_B);

        if(matchToC(fragment))
            genes.add(HLA_C);

        return new NucleotideFragment(fragment.getId(), genes, fragment.getNucleotideLoci(), fragment.getNucleotideQuality(), fragment.getNucleotides());
    }

    private final boolean matchToA(final NucleotideFragment fragment)
    {
        if(fragment.getGenes().contains(HLA_A))
            return true;

        if(fragment.getGenes().contains(HLA_B) && fragment.maxLoci() < 3 * mAbMinBoundary)
            return true;

        if(!fragment.getGenes().contains(HLA_C) && fragment.maxLoci() < 3 * mAcMinBoundary)
            return true;

        return false;
    }

    private final boolean matchToB(final NucleotideFragment fragment)
    {
        if(fragment.getGenes().contains(HLA_B))
            return true;

        if(fragment.getGenes().contains(HLA_A) && fragment.maxLoci() < 3 * mAbMinBoundary)
            return true;

        if(!fragment.getGenes().contains(HLA_C) && fragment.maxLoci() < 3 * mBcMinBoundary)
            return true;

        return false;
    }

    private final boolean matchToC(final NucleotideFragment fragment)
    {
        if(fragment.getGenes().contains(HLA_C))
            return true;

        if(fragment.getGenes().contains(HLA_A) && fragment.maxLoci() < 3 * mAcMinBoundary)
            return true;

        if(!fragment.getGenes().contains(HLA_B) && fragment.maxLoci() < 3 * mBcMinBoundary)
            return true;

        return false;
    }

}
