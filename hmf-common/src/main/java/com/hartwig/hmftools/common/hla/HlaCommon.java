package com.hartwig.hmftools.common.hla;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.BasePosition;

public final class HlaCommon
{
    public static final String HLA_CHROMOSOME_V37 = "6";
    public static final String HLA_CHROMOSOME_V38 = RefGenomeFunctions.enforceChrPrefix(HLA_CHROMOSOME_V37);

    public static final List<String> HLA_GENES = Lists.newArrayList("HLA-A","HLA-B","HLA-C","HLA-DQA1","HLA-DQB1","HLA-DRB1");

    public static final List<GeneData> HLA_GENE_DATA = Lists.newArrayList();

    public static String hlaChromosome(final RefGenomeVersion version) { return version == V37 ? HLA_CHROMOSOME_V37 : HLA_CHROMOSOME_V38; }

    public static void populateGeneData(final List<GeneData> geneDataList)
    {
        HLA_GENE_DATA.addAll(geneDataList.stream()
                .filter(x -> HLA_GENES.contains(x.GeneName)).collect(Collectors.toList()));
    }

    public static boolean containsPosition(final GenomePosition position)
    {
        return containsPosition(position.chromosome(), position.position());
    }

    public static boolean containsPosition(final BasePosition position)
    {
        return containsPosition(position.Chromosome, position.Position);
    }

    public static boolean containsPosition(final String chromosome, final int position)
    {
        if(!chromosome.equals(HLA_CHROMOSOME_V37) && !chromosome.equals(HLA_CHROMOSOME_V38))
            return false;

        return HLA_GENE_DATA.stream()
                .anyMatch(x -> chromosome.equals(x.Chromosome) && positionWithin(position, x.GeneStart, x.GeneEnd));
    }

    public static boolean overlaps(final String chromosome, final int posStart, final int posEnd)
    {
        if(!chromosome.equals(HLA_CHROMOSOME_V37) && !chromosome.equals(HLA_CHROMOSOME_V38))
            return false;

        return HLA_GENE_DATA.stream()
                .anyMatch(x -> chromosome.equals(x.Chromosome) && positionsOverlap(posStart, posEnd, x.GeneStart, x.GeneEnd));
    }
}
