package com.hartwig.hmftools.linx.visualiser.data;

import static com.hartwig.hmftools.linx.visualiser.data.VisExons.sortedDownstreamExons;
import static com.hartwig.hmftools.linx.visualiser.data.VisExons.sortedUpstreamExons;
import static com.hartwig.hmftools.linx.visualiser.file.VisGeneAnnotationType.EXON_LOST;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.linx.visualiser.file.VisFusion;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExon;

public class DisruptedExons
{
    public static List<GenomeRegion> disruptedGeneRegions(final List<VisFusion> fusions, final List<VisGeneExon> exons)
    {
        List<GenomeRegion> result = exons.stream().filter(x -> x.AnnotationType == EXON_LOST)
                .map(x -> GenomeRegions.create(x.Chromosome, x.ExonStart, x.ExonEnd))
                .collect(Collectors.toList());

        for (VisFusion fusion : fusions)
        {
            result.addAll(disruptedGeneRegions(fusion, exons));
        }

        return result;
    }

    public static List<GenomeRegion> disruptedGeneRegions(final VisFusion fusion, final List<VisGeneExon> exons)
    {
        final List<VisGeneExon> upStreamExons = sortedUpstreamExons(fusion, exons).stream()
                .filter(x -> x.ExonRank >= fusion.FusedExonUp).collect(Collectors .toList());

        final List<VisGeneExon> downStreamExons = sortedDownstreamExons(fusion, exons).stream()
                .filter(x -> x.ExonRank <= fusion.FusedExonDown).collect(Collectors.toList());

        if (upStreamExons.isEmpty() || downStreamExons.isEmpty())
        {
            return Collections.emptyList();
        }

        final VisGeneExon finalIncludedUpExon = upStreamExons.get(0);
        final VisGeneExon finalExcludedUpExon = upStreamExons.get(upStreamExons.size() - 1);
        final GenomeRegion upGeneRegion = upGeneExcludedRegion(fusion, finalIncludedUpExon, finalExcludedUpExon);

        final VisGeneExon firstExcludedDownExon = downStreamExons.get(0);
        final VisGeneExon firstIncludedDownExon = downStreamExons.get(downStreamExons.size() - 1);
        final GenomeRegion downGeneRegion = downGeneExcludedRegion(fusion, firstExcludedDownExon, firstIncludedDownExon);

        return Lists.newArrayList(upGeneRegion, downGeneRegion);
    }

    private static Gene downGeneExcludedRegion(final VisFusion fusion, final VisGeneExon firstExcludedDownExon,
            final VisGeneExon firstIncludedDoneExon)
    {
        return fusion.StrandDown < 0 ?
                ImmutableGene.builder()
                        .type(EXON_LOST)
                        .chromosome(firstExcludedDownExon.Chromosome)
                        .end(firstExcludedDownExon.ExonEnd)
                        .start(Math.min(fusion.PosDown, firstIncludedDoneExon.ExonEnd))
                        .strand(fusion.StrandUp)
                        .name(fusion.GeneNameDown)
                        .transcript(fusion.TranscriptDown)
                        .namePosition(0)
                        .build() :
                ImmutableGene.builder()
                        .type(EXON_LOST)
                        .chromosome(firstExcludedDownExon.Chromosome)
                        .start(firstExcludedDownExon.ExonStart)
                        .end(Math.max(fusion.PosDown, firstIncludedDoneExon.ExonStart))
                        .strand(fusion.StrandUp)
                        .name(fusion.GeneNameDown)
                        .transcript(fusion.TranscriptDown)
                        .namePosition(0)
                        .build();
    }

    private static Gene upGeneExcludedRegion(final VisFusion fusion, final VisGeneExon finalIncludedExon,
            final VisGeneExon finalExcludedUpExon)
    {
        return fusion.StrandUp < 0 ?
                ImmutableGene.builder()
                        .type(EXON_LOST)
                        .chromosome(finalExcludedUpExon.Chromosome)
                        .start(finalExcludedUpExon.ExonStart)
                        .end(Math.max(fusion.PosUp, finalIncludedExon.ExonStart))
                        .strand(fusion.StrandDown)
                        .name(fusion.GeneNameUp)
                        .transcript(fusion.TranscriptUp)
                        .namePosition(0)
                        .build() :
                ImmutableGene.builder()
                        .type(EXON_LOST)
                        .chromosome(finalExcludedUpExon.Chromosome)
                        .start(Math.min(fusion.PosUp, finalIncludedExon.ExonEnd))
                        .end(finalExcludedUpExon.ExonEnd)
                        .strand(fusion.StrandDown)
                        .name(fusion.GeneNameUp)
                        .transcript(fusion.TranscriptUp)
                        .namePosition(0)
                        .build();
    }

}
