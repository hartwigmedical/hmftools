package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.lilac.hla.HlaGene_.HLA_A;

import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.hla.HlaGene_;

public enum GeneSelector
{
    // TODO: implement "all" selector

    MHC_CLASS_1((final HlaGene_ gene) -> gene.mhcClass() == MhcClass_.CLASS_1 && !gene.isPseudo()),
    HLA_DQB1(HlaGene_.HLA_DQB1),
    HLA_DPA1(HlaGene_.HLA_DPA1),
    HLA_DPB1(HlaGene_.HLA_DPB1),
    HLA_DQA1(HlaGene_.HLA_DQA1),
    HLA_DRB1(HlaGene_.HLA_DRB1),
    HLA_DRB(List.of(HlaGene_.HLA_DRB1, HlaGene_.HLA_DRB3, HlaGene_.HLA_DRB4, HlaGene_.HLA_DRB5));

    private final LinkedHashSet<HlaGene_> mGenes_;

    GeneSelector(final HlaGene_ gene)
    {
        mGenes_ = Sets.newLinkedHashSet(List.of(gene));
    }

    GeneSelector(final Iterable<HlaGene_> genes) { mGenes_ = Sets.newLinkedHashSet(genes); }

    GeneSelector(final Predicate<HlaGene_> geneFilter)
    {
        mGenes_ = Arrays.stream(HlaGene_.values()).filter(x -> !x.isDebug() && geneFilter.test(x)).collect(Collectors.toCollection(Sets::newLinkedHashSet));
    }

    public LinkedHashSet<HlaGene_> genes_() { return mGenes_; }

    public boolean contains(final HlaGene_ gene)
    {
        return mGenes_.contains(gene);
    }

    public boolean coversMhcClass1()
    {
        return contains(HLA_A);
    }
}
