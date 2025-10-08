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
    MHC_CLASS_1((final HlaGene gene) -> gene.mhcClass() == GeneClass.MHC_CLASS_1 && !gene.isPseudo(), "6"),
    HLA_DQB1(HlaGene.HLA_DQB1, "6"),
    HLA_DPA1(HlaGene.HLA_DPA1, "6"),
    HLA_DPB1(HlaGene.HLA_DPB1, "6"),
    HLA_DQA1(HlaGene.HLA_DQA1, "6"),
    HLA_DRB1(HlaGene.HLA_DRB1, "6"),
    HLA_DRB(List.of(HlaGene.HLA_DRB1, HlaGene.HLA_DRB3, HlaGene.HLA_DRB4, HlaGene.HLA_DRB5), "6"),

    DPYD(HlaGene.DPYD, "1");

    private final LinkedHashSet<HlaGene> mGenes;

    public final String Chromosome;

    GeneSelector(final HlaGene gene, final String chromosome)
    {
        mGenes = Sets.newLinkedHashSet(List.of(gene));
        Chromosome = chromosome;
    }

    GeneSelector(final Predicate<HlaGene> geneFilter, final String chromosome)
    {
        mGenes = Arrays.stream(HlaGene.values())
                .filter(x -> !x.isDebug() && geneFilter.test(x))
                .collect(Collectors.toCollection(Sets::newLinkedHashSet));
        Chromosome = chromosome;
    }

    GeneSelector(final Iterable<HlaGene> genes, final String chromosome)
    {
        mGenes = Sets.newLinkedHashSet(genes);
        Chromosome = chromosome;
    }

    public LinkedHashSet<HlaGene> genes()
    {
        return mGenes;
    }

    public boolean contains(final HlaGene gene)
    {
        return mGenes.contains(gene);
    }

    public boolean coversMhcClass1()
    {
        return contains(HLA_A);
    }
}
