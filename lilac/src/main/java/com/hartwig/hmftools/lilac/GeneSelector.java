package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_A;

import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.hla.HlaGene;

public enum GeneSelector
{
    MHC_CLASS_1(
            (final HlaGene gene) -> gene.mhcClass() == MhcClass.CLASS_1 && !gene.isPseudo(),
            (final HlaGene gene) -> gene.mhcClass() == MhcClass.CLASS_1 && gene.isPseudo()),

    HLA_DQB1(HlaGene.HLA_DQB1),
    HLA_DPA1(HlaGene.HLA_DPA1),
    HLA_DPB1(HlaGene.HLA_DPB1),
    HLA_DQA1(HlaGene.HLA_DQA1),
    HLA_DRB1(HlaGene.HLA_DRB1),
    HLA_DRB(List.of(HlaGene.HLA_DRB1, HlaGene.HLA_DRB3, HlaGene.HLA_DRB4, HlaGene.HLA_DRB5));

    private final LinkedHashSet<HlaGene> mGenes;
    private final LinkedHashSet<HlaGene> mAllGenes;

    GeneSelector(final HlaGene gene)
    {
        mGenes = Sets.newLinkedHashSet(List.of(gene));
        mAllGenes = mGenes;
    }

    GeneSelector(final Predicate<HlaGene> geneFilter, final Predicate<HlaGene> pseudoGeneFilter)
    {
        mGenes = Arrays.stream(HlaGene.values())
                .filter(x -> !x.isDebug() && geneFilter.test(x))
                .collect(Collectors.toCollection(Sets::newLinkedHashSet));

        LinkedHashSet<HlaGene> pseudoGenes = Arrays.stream(HlaGene.values())
                .filter(x -> !x.isDebug() && pseudoGeneFilter.test(x))
                .collect(Collectors.toCollection(Sets::newLinkedHashSet));

        mAllGenes = Sets.newLinkedHashSet();
        mAllGenes.addAll(mGenes);
        mAllGenes.addAll(pseudoGenes);
    }

    GeneSelector(final Iterable<HlaGene> genes)
    {
        mGenes = Sets.newLinkedHashSet(genes);
        mAllGenes = mGenes;
    }

    public LinkedHashSet<HlaGene> genes() { return mGenes; }
    public LinkedHashSet<HlaGene> allGenes() { return mAllGenes; }

    public boolean contains(final HlaGene gene)
    {
        return mGenes.contains(gene);
    }

    public boolean coversMhcClass1()
    {
        return contains(HLA_A);
    }
}
