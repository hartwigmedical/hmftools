package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.HLA_A;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_B;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_C;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;

import org.apache.commons.math3.util.Pair;

public class NucleotideGeneEnrichment
{
    private final int mAbMinBoundary;
    private final int mAcMinBoundary;
    private final int mBcMinBoundary;

    public NucleotideGeneEnrichment(final List<Integer> aBoundaries, final List<Integer> bBoundaries, final List<Integer> cBoundaries)
    {
        Set<Integer> abUnique = Sets.newHashSet();
        Set<Integer> acUnique = Sets.newHashSet();
        Set<Integer> bcUnique = Sets.newHashSet();

        aBoundaries.stream().filter(x -> !bBoundaries.contains(x)).forEach(x -> abUnique.add(x));
        bBoundaries.stream().filter(x -> !aBoundaries.contains(x)).forEach(x -> abUnique.add(x));
        mAbMinBoundary = abUnique.stream().mapToInt(x -> x).min().orElse(0);

        aBoundaries.stream().filter(x -> !cBoundaries.contains(x)).forEach(x -> acUnique.add(x));
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

    public List<NucleotideFragment> enrich(final List<NucleotideFragment> fragments)
    {
        return fragments.stream().map(x -> enrich(x)).collect(Collectors.toList());
    }

    public final NucleotideFragment enrich(final NucleotideFragment fragment)
    {
        if(fragment.containsIndel())
            return fragment;

        Set<String> genes = Sets.newHashSet();

        if(matchToGene(fragment, HLA_A, new Pair(HLA_B, mAbMinBoundary), new Pair(HLA_C, mAcMinBoundary)))
            genes.add(HLA_A);

        if(matchToGene(fragment, HLA_B, new Pair(HLA_A, mAbMinBoundary), new Pair(HLA_C, mBcMinBoundary)))
            genes.add(HLA_B);

        if(matchToGene(fragment, HLA_C, new Pair(HLA_A, mAcMinBoundary), new Pair(HLA_B, mBcMinBoundary)))
            genes.add(HLA_C);

        return new NucleotideFragment(fragment.getId(), genes, fragment.getNucleotideLoci(), fragment.getNucleotideQuality(), fragment.getNucleotides());
    }

    private final boolean matchToGene(
            final NucleotideFragment fragment, final String gene,
            final Pair<String,Integer> otherGene1, final Pair<String,Integer> otherGene2)
    {
        if(fragment.containsGene(gene))
            return true;

        if(fragment.containsGene(otherGene1.getFirst()) && fragment.maxLoci() < 3 * otherGene1.getSecond())
            return true;

        if(fragment.containsGene(otherGene2.getFirst()) && fragment.maxLoci() < 3 * otherGene2.getSecond())
            return true;

        return false;
    }
}
