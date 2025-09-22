package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;
import static com.hartwig.hmftools.lilac.hla.HlaGene_.HLA_A;
import static com.hartwig.hmftools.lilac.hla.HlaGene_.HLA_DRB3;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.OptionalInt;
import java.util.Set;

import com.beust.jcommander.internal.Lists;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.GeneSelector;
import com.hartwig.hmftools.lilac.hla.HlaGene_;

public class NucleotideGeneEnrichment
{
    private final List<HlaGene_> mGenes;
    private final HashMap<Set<HlaGene_>, OptionalInt> mMinUniqueProteinExonBoundaries;

    private NucleotideGeneEnrichment(final List<HlaGene_> genes, final HashMap<Set<HlaGene_>, OptionalInt> minUniqueProteinExonBoundaries)
    {
        mGenes = genes;
        mMinUniqueProteinExonBoundaries = minUniqueProteinExonBoundaries;
    }

    public static NucleotideGeneEnrichment create(final Map<HlaGene_, List<Integer>> geneBoundaries_)
    {
        // determine the minimum unique exon boundary for each pair
        List<HlaGene_> genes = null;
        if(geneBoundaries_.containsKey(HLA_A))
        {
            genes = Lists.newArrayList(GeneSelector.MHC_CLASS_1.genes_());
        }
        else if(geneBoundaries_.containsKey(HLA_DRB3))
        {
            genes = Lists.newArrayList(GeneSelector.HLA_DRB.genes_());
        }
        else
        {
            return null;
        }

        HashMap<Set<HlaGene_>, OptionalInt> minUniqueProteinExonBoundaries = Maps.newHashMap();
        for(int i = 0; i < genes.size() - 1; i++)
        {
            HlaGene_ gene1 = genes.get(i);
            for(int j = i + 1; j < genes.size(); j++)
            {
                HlaGene_ gene2 = genes.get(j);
                Set<HlaGene_> key = Sets.newHashSet(gene1, gene2);
                OptionalInt minUniqueExonBoundary = getMinUniqueBoundary(geneBoundaries_.get(gene1), geneBoundaries_.get(gene2));
                minUniqueProteinExonBoundaries.put(key, minUniqueExonBoundary);
            }
        }

        return new NucleotideGeneEnrichment(genes, minUniqueProteinExonBoundaries);
    }

    private static OptionalInt getMinUniqueBoundary(final List<Integer> boundariesGene1, final List<Integer> boundariesGene2)
    {
        Set<Integer> uniqueBoundaries = Sets.newHashSet();

        boundariesGene1.stream().filter(x -> !boundariesGene2.contains(x)).forEach(uniqueBoundaries::add);
        boundariesGene2.stream().filter(x -> !boundariesGene1.contains(x)).forEach(uniqueBoundaries::add);

        if(uniqueBoundaries.isEmpty())
            return OptionalInt.empty();

        return OptionalInt.of(uniqueBoundaries.stream().mapToInt(x -> x).min().orElse(0));
    }

    public NavigableMap<Integer, Integer> getFilters(final HlaGene_ gene)
    {
        NavigableMap<Integer, Integer> filters = Maps.newTreeMap();
        for(HlaGene_ otherGene : mGenes)
        {
            if(gene == otherGene)
                continue;

            OptionalInt filter = getFilter(gene, otherGene);
            if(filter.isEmpty())
                continue;

            filters.merge(filter.getAsInt(), 1, Integer::sum);
        }

        return filters;
    }

    public OptionalInt getFilter(final HlaGene_ gene1, final HlaGene_ gene2)
    {
        Set<HlaGene_> key = Sets.newHashSet(gene1, gene2);
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
        for(HlaGene_ gene : mGenes)
        {
            if(!considerAddingGene(fragment, gene, maxFragmentNucleotideLocus))
                continue;

            boolean addGene = false;
            for(HlaGene_ otherGene : mGenes)
            {
                if(gene == otherGene)
                    continue;

                OptionalInt minUniqueProteinExonBoundary = getFilter(gene, otherGene);
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

    private static boolean considerAddingGene(final Fragment fragment, final HlaGene_ gene, int maxFragmentNucleotideLocus)
    {
        if(fragment.containsGene(gene))
            return false;

        // the max supported read locus cannot be past the end of the gene's range
        return maxFragmentNucleotideLocus <= GENE_CACHE.NucleotideLengths.get(gene);
    }

    private static boolean checkAddAdditionalGene(
            final Fragment fragment, int maxFragmentNucleotideLocus, final HlaGene_ otherGene, final OptionalInt geneComboUniqueAminoAcidBoundary)
    {
        if(!fragment.containsGene(otherGene))
            return false;

        return geneComboUniqueAminoAcidBoundary.isEmpty() || maxFragmentNucleotideLocus < geneComboUniqueAminoAcidBoundary.getAsInt() * 3;
    }
}
