package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_A;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_B;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_C;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.hla.HlaGene;

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

    private static int getMinUniqueBoundary(final List<Integer> boundariesGene1, final List<Integer> boundariesGene2)
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

    public void checkAddAdditionalGenes(final List<Fragment> fragments)
    {
        fragments.forEach(x -> checkAdditionalGenes(x));
    }

    @VisibleForTesting
    public void checkAdditionalGenes(final Fragment fragment)
    {
        if(fragment.containsIndel())
            return;

        // logic: for any gene which isn't yet associated with the fragment, test if it could be by check its max nucleotide locus
        // versus the first unique nucelotide for the gene pair
        //
        // example: a fragment isn't associated with gene A, is with A, and the fragments max base is within the unique base of A and B,
        // so cannot it be distinguished between them - then add it to A as well
        int maxFragmentNucleotideLocus = fragment.maxNucleotideLocus();
        if(considerAddingGene(fragment, HLA_A, maxFragmentNucleotideLocus))
        {
            if(checkAddAdditionalGene(fragment, maxFragmentNucleotideLocus, HLA_B, mAbMinUniqueProteinExonBoundary)
            || checkAddAdditionalGene(fragment, maxFragmentNucleotideLocus, HLA_C, mAcMinUniqueProteinExonBoundary))
            {
                fragment.addGene(HLA_A);
            }
        }

        if(considerAddingGene(fragment, HLA_B, maxFragmentNucleotideLocus))
        {
            // test: HLA_B, mAbMinUniqueProteinExonBoundary, HLA_C, mAcMinUniqueProteinExonBoundary))
            if(checkAddAdditionalGene(fragment, maxFragmentNucleotideLocus, HLA_A, mAbMinUniqueProteinExonBoundary)
            || checkAddAdditionalGene(fragment, maxFragmentNucleotideLocus, HLA_C, mBcMinUniqueProteinExonBoundary))
            {
                fragment.addGene(HLA_B);
            }
        }

        if(considerAddingGene(fragment, HLA_C, maxFragmentNucleotideLocus))
        {
            // test: HLA_B, mAbMinUniqueProteinExonBoundary, HLA_C, mAcMinUniqueProteinExonBoundary))
            if(checkAddAdditionalGene(fragment, maxFragmentNucleotideLocus, HLA_A, mAcMinUniqueProteinExonBoundary)
            || checkAddAdditionalGene(fragment, maxFragmentNucleotideLocus, HLA_B, mBcMinUniqueProteinExonBoundary))
            {
                fragment.addGene(HLA_C);
            }
        }
    }

    private static boolean considerAddingGene(final Fragment fragment, final HlaGene gene, int maxFragmentNucleotideLocus)
    {
        if(fragment.containsGene(gene))
            return false;

        // the max supported read locus cannot be past the end of the gene's range
        return maxFragmentNucleotideLocus <= GENE_CACHE.NucleotideLengths.get(gene);
    }

    private boolean checkAddAdditionalGene(
            final Fragment fragment, int maxFragmentNucleotideLocus, final HlaGene otherGene, int geneComboUniqueAminoAcidBoundary)
    {
        if(!fragment.containsGene(otherGene))
            return false;

        return maxFragmentNucleotideLocus < geneComboUniqueAminoAcidBoundary * 3;
    }
}
