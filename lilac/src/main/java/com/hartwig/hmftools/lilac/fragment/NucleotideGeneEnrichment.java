package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_A;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_B;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_C;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.MhcClass_;
import com.hartwig.hmftools.lilac.hla.HlaGene;

import org.apache.commons.lang3.NotImplementedException;

public class NucleotideGeneEnrichment
{
    private final int mAbMinUniqueProteinExonBoundary_;
    private final int mAcMinUniqueProteinExonBoundary_;
    private final int mBcMinUniqueProteinExonBoundary_;

    public NucleotideGeneEnrichment(final Map<HlaGene, List<Integer>> geneBoundaries_)
    {
        // determine the minimum unique exon boundary for each pair
        mAbMinUniqueProteinExonBoundary_ = getMinUniqueBoundary(geneBoundaries_.get(HLA_A), geneBoundaries_.get(HLA_B));
        mAcMinUniqueProteinExonBoundary_ = getMinUniqueBoundary(geneBoundaries_.get(HLA_A), geneBoundaries_.get(HLA_C));
        mBcMinUniqueProteinExonBoundary_ = getMinUniqueBoundary(geneBoundaries_.get(HLA_B), geneBoundaries_.get(HLA_C));
    }

    private static int getMinUniqueBoundary(final List<Integer> boundariesGene1, final List<Integer> boundariesGene2)
    {
        Set<Integer> uniqueBoundaries = Sets.newHashSet();

        boundariesGene1.stream().filter(x -> !boundariesGene2.contains(x)).forEach(uniqueBoundaries::add);
        boundariesGene2.stream().filter(x -> !boundariesGene1.contains(x)).forEach(uniqueBoundaries::add);

        return uniqueBoundaries.stream().mapToInt(x -> x).min().orElse(0);
    }

    public int getAFilterB_() { return mAbMinUniqueProteinExonBoundary_; }
    public int getAFilterC_() { return mAcMinUniqueProteinExonBoundary_; }
    public int getBFilterA_() { return mAbMinUniqueProteinExonBoundary_; }
    public int getBFilterC_() { return mBcMinUniqueProteinExonBoundary_; }
    public int getCFilterA_() { return mAcMinUniqueProteinExonBoundary_; }
    public int getCFilterB_() { return mBcMinUniqueProteinExonBoundary_; }

    public void checkAddAdditionalGenes(final Iterable<Fragment> fragments)
    {
        fragments.forEach(this::checkAdditionalGenes);
    }

    @VisibleForTesting
    public void checkAdditionalGenes(final Fragment fragment)
    {
        if(fragment.containsIndel())
            return;

        if(fragment.readGene().mhcClass() != MhcClass_.CLASS_1)
            return;

        // logic: for any gene which isn't yet associated with the fragment, test if it could be by check its max nucleotide locus
        // versus the first unique nucelotide for the gene pair
        //
        // example: a fragment isn't associated with gene A, is with A, and the fragments max base is within the unique base of A and B,
        // so cannot it be distinguished between them - then add it to A as well
        int maxFragmentNucleotideLocus = fragment.maxNucleotideLocus();
        if(considerAddingGene(fragment, HLA_A, maxFragmentNucleotideLocus))
        {
            if(checkAddAdditionalGene(fragment, maxFragmentNucleotideLocus, HLA_B, mAbMinUniqueProteinExonBoundary_)
                    || checkAddAdditionalGene(fragment, maxFragmentNucleotideLocus, HLA_C, mAcMinUniqueProteinExonBoundary_))
            {
                fragment.addGene(HLA_A);
            }
        }

        if(considerAddingGene(fragment, HLA_B, maxFragmentNucleotideLocus))
        {
            // test: HLA_B, mAbMinUniqueProteinExonBoundary, HLA_C, mAcMinUniqueProteinExonBoundary))
            if(checkAddAdditionalGene(fragment, maxFragmentNucleotideLocus, HLA_A, mAbMinUniqueProteinExonBoundary_)
                    || checkAddAdditionalGene(fragment, maxFragmentNucleotideLocus, HLA_C, mBcMinUniqueProteinExonBoundary_))
            {
                fragment.addGene(HLA_B);
            }
        }

        if(considerAddingGene(fragment, HLA_C, maxFragmentNucleotideLocus))
        {
            // test: HLA_B, mAbMinUniqueProteinExonBoundary, HLA_C, mAcMinUniqueProteinExonBoundary))
            if(checkAddAdditionalGene(fragment, maxFragmentNucleotideLocus, HLA_A, mAcMinUniqueProteinExonBoundary_)
                    || checkAddAdditionalGene(fragment, maxFragmentNucleotideLocus, HLA_B, mBcMinUniqueProteinExonBoundary_))
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

    private static boolean checkAddAdditionalGene(
            final Fragment fragment, int maxFragmentNucleotideLocus, final HlaGene otherGene, int geneComboUniqueAminoAcidBoundary)
    {
        if(!fragment.containsGene(otherGene))
            return false;

        return maxFragmentNucleotideLocus < geneComboUniqueAminoAcidBoundary * 3;
    }
}
