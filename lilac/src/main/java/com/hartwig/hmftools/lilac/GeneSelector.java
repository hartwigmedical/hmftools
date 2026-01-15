package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_A;

import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
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
    HLA_DRB1(HlaGene.HLA_DRB1, false),
    HLA_DRB(List.of(HlaGene.HLA_DRB1, HlaGene.HLA_DRB3, HlaGene.HLA_DRB4, HlaGene.HLA_DRB5), false);

    public static final String ALL_GENES = "ALL";

    public final boolean IncludedInAll;

    private final LinkedHashSet<HlaGene> mGenes;
    private final LinkedHashSet<HlaGene> mAllGenes;

    GeneSelector(final HlaGene gene)
    {
        this(gene, true);
    }

    GeneSelector(final HlaGene gene, boolean includedInAll)
    {
        IncludedInAll = includedInAll;
        mGenes = Sets.newLinkedHashSet(List.of(gene));
        mAllGenes = mGenes;
    }

    GeneSelector(final Predicate<HlaGene> geneFilter, final Predicate<HlaGene> pseudoGeneFilter)
    {
        IncludedInAll = true;
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
        this(genes, true);
    }

    GeneSelector(final Iterable<HlaGene> genes, boolean includedInAll)
    {
        IncludedInAll = includedInAll;
        mGenes = Sets.newLinkedHashSet(genes);
        mAllGenes = mGenes;
    }

    public LinkedHashSet<HlaGene> genes() { return mGenes; }
    public int geneCount() { return mGenes.size(); }
    public LinkedHashSet<HlaGene> allGenes() { return mAllGenes; }

    public boolean contains(final HlaGene gene)
    {
        return mGenes.contains(gene);
    }

    public boolean coversMhcClass1()
    {
        return contains(HLA_A);
    }

    public static Set<GeneSelector> parseArg(final String arg)
    {
        Set<GeneSelector> genes = Sets.newHashSet();
        String[] geneGroups = arg.split(",");
        for(String geneGroup : geneGroups)
        {
            if(geneGroup.equalsIgnoreCase(ALL_GENES))
            {
                for(GeneSelector geneSelector : values())
                {
                    if(geneSelector.IncludedInAll)
                        genes.add(geneSelector);
                }

                continue;
            }

            genes.add(valueOf(geneGroup));
        }

        return genes;
    }
}
