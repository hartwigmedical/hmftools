package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.HLA_A;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_B;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_C;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;

public class NucleotideGeneEnrichment
{
    private final int mAbMinUniqueProteinExonBoundary;
    private final int mAcMinUniqueProteinExonBoundary;
    private final int mBcMinUniqueProteinExonBoundary;

    public NucleotideGeneEnrichment(final List<Integer> aBoundaries, final List<Integer> bBoundaries, final List<Integer> cBoundaries)
    {
        // determine the minimum unique exon boundary for each pair
        mAbMinUniqueProteinExonBoundary = getMinUniqueBoundary(aBoundaries, bBoundaries);
        mAcMinUniqueProteinExonBoundary = getMinUniqueBoundary(aBoundaries, cBoundaries);
        mBcMinUniqueProteinExonBoundary = getMinUniqueBoundary(bBoundaries, cBoundaries);
    }

    private static int getMinUniqueBoundary(List<Integer> boundariesGene1, List<Integer> boundariesGene2)
    {
        Set<Integer> uniqueBoundaries = Sets.newHashSet();

        boundariesGene1.stream().filter(x -> !boundariesGene2.contains(x)).forEach(x -> uniqueBoundaries.add(x));
        boundariesGene2.stream().filter(x -> !boundariesGene1.contains(x)).forEach(x -> uniqueBoundaries.add(x));

        return uniqueBoundaries.stream().mapToInt(x -> x).min().orElse(0);
    }

    public final int getAFilterB() { return mAbMinUniqueProteinExonBoundary; }
    public final int getAFilterC() { return mAcMinUniqueProteinExonBoundary; }
    public final int getBFilterA() { return mAbMinUniqueProteinExonBoundary; }
    public final int getBFilterC() { return mBcMinUniqueProteinExonBoundary; }
    public final int getCFilterA() { return mAcMinUniqueProteinExonBoundary; }
    public final int getCFilterB() { return mBcMinUniqueProteinExonBoundary; }

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
        if(!fragment.containsGene(HLA_A) && matchToGenes(fragment, Map.of(HLA_B, mAbMinUniqueProteinExonBoundary, HLA_C, mAcMinUniqueProteinExonBoundary)))
            fragment.addGene(HLA_A);

        if(!fragment.containsGene(HLA_B) && matchToGenes(fragment, Map.of(HLA_A, mAbMinUniqueProteinExonBoundary, HLA_C, mBcMinUniqueProteinExonBoundary)))
            fragment.addGene(HLA_B);

        if(!fragment.containsGene(HLA_C) && matchToGenes(fragment, Map.of(HLA_A, mAcMinUniqueProteinExonBoundary, HLA_B, mBcMinUniqueProteinExonBoundary)))
            fragment.addGene(HLA_C);

        return fragment;
    }

    private boolean matchToGenes(final Fragment fragment, final Map<String,Integer> minUniqueProteinExonBoundariesPerGene)
    {
        for(String geneName : minUniqueProteinExonBoundariesPerGene.keySet())
        {
            int minUniqueNucExonBoundary = minUniqueProteinExonBoundariesPerGene.get(geneName) * 3;

            if(fragment.containsGene(geneName) && fragment.maxNucleotideLocus() < minUniqueNucExonBoundary)
                return true;
        }

        return false;
    }

    // TODO: Rewrite matchToGenes() to be compatible with class GeneCache
    // private boolean matchToGenes(
    //         final Fragment fragment,
    //         final Map<String, List<Integer>> proteinExonBoundariesList)
    // {
    //     for(String geneName : proteinExonBoundariesList.keySet())
    //     {
    //         List<Integer> proteinExonBoundaries = proteinExonBoundariesList.get(geneName);
    //         // int minUniqueProteinExonBoundary = getMinUniqueBoundary(proteinExonBoundaries);
    //     }
    // }
}
