package com.hartwig.hmftools.linx.visualiser.circos;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.linx.visualiser.data.FusedExon;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableFusedExon;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableProteinDomain;
import com.hartwig.hmftools.linx.visualiser.data.ProteinDomain;

import org.jetbrains.annotations.NotNull;

class ScaleIntrons
{
    private static final int INTRON_LENGTH = 10;

    private final List<GenomeRegion> introns;

    ScaleIntrons(@NotNull final List<GenomeRegion> introns)
    {
        this.introns = introns;
    }

    @NotNull
    public static List<GenomeRegion> introns(@NotNull final List<? extends GenomeRegion> exons)
    {
        final List<GenomeRegion> result = Lists.newArrayList();
        for (int i = 0; i < exons.size() - 1; i++)
        {
            GenomeRegion current = exons.get(i);
            GenomeRegion next = exons.get(i + 1);

            if (next.start() > current.end())
            {
                result.add(GenomeRegions.create(current.chromosome(), current.end() + 1, next.start() - 1));
            }
        }

        return result;
    }

    @NotNull
    public List<FusedExon> scaleIntronsInExons(@NotNull final List<FusedExon> exons)
    {
        final List<FusedExon> result = Lists.newArrayList();
        for (FusedExon exon : exons)
        {
            result.add(ImmutableFusedExon.builder().from(exon)
                    .start(exon.start() + offset(exon.start()))
                    .end(exon.end() + offset(exon.end()))
                    .geneStart(exon.geneStart() + offset(exon.geneStart()))
                    .geneEnd(exon.geneEnd() + offset(exon.geneEnd()))
                    .build());
        }

        return result;
    }

    @NotNull
    public List<ProteinDomain> scaleIntronsInProteinDomains(@NotNull final List<ProteinDomain> domains)
    {
        final List<ProteinDomain> result = Lists.newArrayList();
        for (ProteinDomain domain : domains)
        {
            result.add(ImmutableProteinDomain.builder().from(domain)
                    .start(domain.start() + offset(domain.start()))
                    .end(domain.end() + offset(domain.end()))
                    .build());
        }

        return result;
    }

    @NotNull
    public GenomeRegion scaleIntronsFromRegion(@NotNull final GenomeRegion region)
    {
        return GenomeRegions.create(region.chromosome(), region.start() + offset(region.start()), region.end() + offset(region.end()));
    }

    private int offset(long position)
    {
        int offset = 0;
        for (GenomeRegion intron : introns)
        {
            if (position > intron.start())
            {
                long distanceFromStart = Math.min(intron.bases(), position - intron.start() + 1);
                double proportion = (1d * distanceFromStart) / intron.bases();
                offset += Math.round(proportion * INTRON_LENGTH) - distanceFromStart;
            }
        }

        return offset;
    }
}
