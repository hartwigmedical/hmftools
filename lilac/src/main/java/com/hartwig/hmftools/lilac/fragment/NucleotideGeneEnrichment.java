package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_A;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_DRB1;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Set;

import com.beust.jcommander.internal.Lists;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.GeneSelector;
import com.hartwig.hmftools.lilac.hla.HlaGene;

public class NucleotideGeneEnrichment
{
    private final List<HlaGene> mGenes;
    private final HashMap<Set<HlaGene>, Integer> mMinUniqueProteinExonBoundaries;

    private NucleotideGeneEnrichment(final List<HlaGene> genes, final HashMap<Set<HlaGene>, Integer> minUniqueProteinExonBoundaries)
    {
        mGenes = genes;
        mMinUniqueProteinExonBoundaries = minUniqueProteinExonBoundaries;
    }

    public static NucleotideGeneEnrichment create(final Map<HlaGene, List<Integer>> geneBoundaries)
    {
        // determine the minimum unique exon boundary for each pair
        List<HlaGene> genes;
        if(geneBoundaries.containsKey(HLA_A))
        {
            genes = Lists.newArrayList(GeneSelector.MHC_CLASS_1.genes());
        }
        else if(geneBoundaries.containsKey(HLA_DRB1))
        {
            genes = Lists.newArrayList(GeneSelector.HLA_DRB.genes());
        }
        else
        {
            return null;
        }

        HashMap<Set<HlaGene>, Integer> minUniqueProteinExonBoundaries = Maps.newHashMap();
        for(int i = 0; i < genes.size() - 1; i++)
        {
            HlaGene gene1 = genes.get(i);
            for(int j = i + 1; j < genes.size(); j++)
            {
                HlaGene gene2 = genes.get(j);
                Set<HlaGene> key = Sets.newHashSet(gene1, gene2);
                Integer minUniqueExonBoundary = getMinUniqueBoundary(geneBoundaries.get(gene1), geneBoundaries.get(gene2));
                minUniqueProteinExonBoundaries.put(key, minUniqueExonBoundary);
            }
        }

        return new NucleotideGeneEnrichment(genes, minUniqueProteinExonBoundaries);
    }

    private static Integer getMinUniqueBoundary(final List<Integer> boundariesGene1, final List<Integer> boundariesGene2)
    {
        Set<Integer> uniqueBoundaries = Sets.newHashSet();

        boundariesGene1.stream().filter(x -> !boundariesGene2.contains(x)).forEach(uniqueBoundaries::add);
        boundariesGene2.stream().filter(x -> !boundariesGene1.contains(x)).forEach(uniqueBoundaries::add);

        if(uniqueBoundaries.isEmpty())
            return null;

        return uniqueBoundaries.stream().mapToInt(x -> x).min().orElse(0);
    }

    public NavigableMap<Integer, Integer> getFilters(final HlaGene gene)
    {
        NavigableMap<Integer, Integer> filters = Maps.newTreeMap();
        for(HlaGene otherGene : mGenes)
        {
            if(gene == otherGene)
                continue;

            Integer filter = getFilter(gene, otherGene);
            if(filter == null)
                continue;

            filters.merge(filter, 1, Integer::sum);
        }

        return filters;
    }

    public Integer getFilter(final HlaGene gene1, final HlaGene gene2)
    {
        Set<HlaGene> key = Sets.newHashSet(gene1, gene2);
        return mMinUniqueProteinExonBoundaries.get(key);
    }

    public void checkAddAdditionalGenes(final Iterable<Fragment> fragments)
    {
        fragments.forEach(this::checkAdditionalGenes);
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
        for(HlaGene gene : mGenes)
        {
            if(!considerAddingGene(fragment, gene, maxFragmentNucleotideLocus))
                continue;

            boolean addGene = false;
            for(HlaGene otherGene : mGenes)
            {
                if(gene == otherGene)
                    continue;

                Integer minUniqueProteinExonBoundary = getFilter(gene, otherGene);
                if(!checkAddAdditionalGene(fragment, maxFragmentNucleotideLocus, otherGene, minUniqueProteinExonBoundary))
                    continue;

                addGene = true;
                break;
            }

            if(!addGene)
                continue;

            fragment.addGene(gene);
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
            final Fragment fragment, int maxFragmentNucleotideLocus, final HlaGene otherGene,
            final Integer geneComboUniqueAminoAcidBoundary)
    {
        if(!fragment.containsGene(otherGene))
            return false;

        return geneComboUniqueAminoAcidBoundary == null || maxFragmentNucleotideLocus < geneComboUniqueAminoAcidBoundary * 3;
    }
}
