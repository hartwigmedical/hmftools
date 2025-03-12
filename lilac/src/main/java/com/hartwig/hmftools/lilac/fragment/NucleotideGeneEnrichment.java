package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.HLA_A;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_B;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_C;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;

import org.apache.commons.math3.util.Pair;

public class NucleotideGeneEnrichment
{
    private final int mAbMinBoundary;
    private final int mAcMinBoundary;
    private final int mBcMinBoundary;

    public NucleotideGeneEnrichment(final List<Integer> aBoundaries, final List<Integer> bBoundaries, final List<Integer> cBoundaries)
    {
        // determine the minimum unique exon boundary for each pair
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

    public List<Fragment> checkAddAdditionalGenes(final List<Fragment> fragments)
    {
        return fragments.stream().map(x -> checkAddAdditionalGenes(x)).collect(Collectors.toList());
    }

    @VisibleForTesting
    public final Fragment checkAddAdditionalGenes(final Fragment fragment)
    {
        if(fragment.containsIndel())
            return fragment;

        // add an extra gene for consideration if the fragment's max loci is within the minimum unique boundary for either pair
        if(!fragment.containsGene(HLA_A) && matchToGene(fragment, Pair.create(HLA_B, mAbMinBoundary), Pair.create(HLA_C, mAcMinBoundary)))
            fragment.addGene(HLA_A);

        if(!fragment.containsGene(HLA_B) && matchToGene(fragment, Pair.create(HLA_A, mAbMinBoundary), Pair.create(HLA_C, mBcMinBoundary)))
            fragment.addGene(HLA_B);

        if(!fragment.containsGene(HLA_C) && matchToGene(fragment, Pair.create(HLA_A, mAcMinBoundary), Pair.create(HLA_B, mBcMinBoundary)))
            fragment.addGene(HLA_C);

        return fragment;
    }

    private boolean matchToGene(
            final Fragment fragment, final Pair<String,Integer> otherGene1, final Pair<String,Integer> otherGene2)
    {
        if(fragment.containsGene(otherGene1.getFirst()) && fragment.maxNucleotideLocus() < 3 * otherGene1.getSecond())
            return true;

        if(fragment.containsGene(otherGene2.getFirst()) && fragment.maxNucleotideLocus() < 3 * otherGene2.getSecond())
            return true;

        return false;
    }
}
