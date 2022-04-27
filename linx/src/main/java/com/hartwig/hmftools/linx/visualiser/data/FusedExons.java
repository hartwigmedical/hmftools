package com.hartwig.hmftools.linx.visualiser.data;

import static com.hartwig.hmftools.linx.visualiser.data.VisExons.sortedDownstreamExons;
import static com.hartwig.hmftools.linx.visualiser.data.VisExons.sortedUpstreamExons;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.linx.visualiser.file.VisFusion;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExon;

import org.jetbrains.annotations.NotNull;

public class FusedExons
{
    private static final String DELIMITER = "\t";

    public static void write(final String fileName, final List<FusedExon> fusedExons) throws IOException
    {
        Files.write(new File(fileName).toPath(), toLines(fusedExons));
    }

    public static List<FusedExon> fusedExons(final VisFusion fusion, final List<VisGeneExon> exons)
    {
        final List<FusedExon> result = Lists.newArrayList();

        final List<VisGeneExon> fusionExons = VisExons.fusionExons(fusion, exons);
        final List<VisGeneExon> upStreamExons = sortedUpstreamExons(fusion, fusionExons);
        final List<VisGeneExon> downStreamExons = sortedDownstreamExons(fusion, fusionExons);

        if (upStreamExons.isEmpty() || downStreamExons.isEmpty())
            return result;

        final VisGeneExon firstUpExon = upStreamExons.get(0);
        final GenomeRegion upGeneRegion = upGeneRegion(fusion, firstUpExon);
        final GenomeRegion convertedUpGeneRegion = convertRegion(fusion.StrandUp, upGeneRegion, upGeneRegion);

        final ImmutableFusedExon.Builder upFusedExonBuilder = ImmutableFusedExon.builder()
                .sampleId(fusion.SampleId)
                .clusterId(fusion.ClusterId)
                .fusion(fusion.name())
                .chromosome(fusion.ChrUp)
                .unadjustedGeneStart(upGeneRegion.start())
                .gene(fusion.GeneNameUp)
                .geneStart(convertedUpGeneRegion.start())
                .geneEnd(convertedUpGeneRegion.end())
                .transcript(fusion.TranscriptUp);

        for (final VisGeneExon exon : upStreamExons)
        {
            final GenomeRegion convertedExon = convertRegion(fusion.StrandUp, upGeneRegion, exon);

            if(upGeneRegion.overlaps(exon))
            {
                final FusedExon fusedExon = upFusedExonBuilder
                        .start(convertedExon.start())
                        .end(convertedExon.end())
                        .rank(exon.ExonRank)
                        .skipped(exon.ExonRank > fusion.FusedExonUp)
                        .build();
                result.add(fusedExon);
            }
        }

        final VisGeneExon finalDownExon = downStreamExons.get(downStreamExons.size() - 1);
        final GenomeRegion downGeneRegion = downGene(fusion, finalDownExon);
        final GenomeRegion convertedDownGeneRegion = convertRegion(fusion.StrandDown, downGeneRegion, downGeneRegion);

        final ImmutableFusedExon.Builder downFusedExonBuilder = ImmutableFusedExon.builder().from(upFusedExonBuilder.build())
                .chromosome(fusion.ChrDown)
                .unadjustedGeneStart(fusion.PosDown)
                .gene(fusion.GeneNameDown)
                .geneStart(convertedDownGeneRegion.start() + convertedUpGeneRegion.end())
                .geneEnd(convertedDownGeneRegion.end() + convertedUpGeneRegion.end())
                .transcript(fusion.TranscriptDown);

        boolean intronicToExonicFusion = fusion.RegionTypeUp.equals("Intronic") && fusion.RegionTypeDown.equals("Exonic");

        for (int i = 0; i < downStreamExons.size(); i++)
        {
            final VisGeneExon exon = downStreamExons.get(i);
            final GenomeRegion convertedExon = convertRegion(fusion.StrandDown, downGeneRegion, exon);

            if(downGeneRegion.overlaps(exon))
            {
                final FusedExon fusedExon = downFusedExonBuilder
                        .start(convertedExon.start() + convertedUpGeneRegion.end())
                        .end(convertedExon.end() + convertedUpGeneRegion.end())
                        .rank(exon.ExonRank)
                        .skipped(exon.ExonRank < fusion.FusedExonDown || (i == 0 && intronicToExonicFusion))
                        .build();
                result.add(fusedExon);
            }
        }

        return result;
    }

    static List<String> toLines(final List<FusedExon> exons)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        exons.stream().map(FusedExons::toString).forEach(lines::add);
        return lines;
    }
    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER).add("sampleId")
                .add("clusterId")
                .add("fusion")
                .add("gene")
                .add("geneStart")
                .add("geneEnd")
                .add("chromosome")
                .add("start")
                .add("end")
                .add("rank")
                .add("skipped")
                .add("transcript")
                .toString();
    }

    @NotNull
    private static String toString(final FusedExon exon)
    {
        return new StringJoiner(DELIMITER)
                .add(exon.sampleId())
                .add(String.valueOf(exon.clusterId()))
                .add(exon.fusion())
                .add(exon.gene())
                .add(String.valueOf(exon.geneStart()))
                .add(String.valueOf(exon.geneEnd()))
                .add(exon.chromosome())
                .add(String.valueOf(exon.start()))
                .add(String.valueOf(exon.end()))
                .add(String.valueOf(exon.rank()))
                .add(String.valueOf(exon.skipped()))
                .add(exon.transcript())
                .toString();
    }

    private static Gene upGeneRegion(final VisFusion fusion, final VisGeneExon firstUpExon)
    {
        return fusion.StrandUp < 0 ?
                ImmutableGene.builder()
                        .type(firstUpExon.AnnotationType)
                        .chromosome(firstUpExon.Chromosome)
                        .start(fusion.PosUp)
                        .end(firstUpExon.ExonEnd)
                        .strand(fusion.StrandUp)
                        .name(fusion.GeneNameDown)
                        .transcript(fusion.TranscriptUp)
                        .namePosition(0)
                        .build() :
                ImmutableGene.builder()
                        .type(firstUpExon.AnnotationType)
                        .chromosome(firstUpExon.Chromosome)
                        .start(firstUpExon.ExonStart)
                        .end(fusion.PosUp)
                        .strand(fusion.StrandUp)
                        .name(fusion.GeneNameDown)
                        .transcript(fusion.TranscriptDown)
                        .namePosition(0)
                        .build();
    }

    static Gene downGene(final VisFusion fusion, final VisGeneExon finalDownGene)
    {
        return fusion.StrandDown < 0 ?
                ImmutableGene.builder()
                        .type(finalDownGene.AnnotationType)
                        .chromosome(finalDownGene.Chromosome)
                        .start(finalDownGene.ExonStart)
                        .end(fusion.PosDown)
                        .strand(fusion.StrandDown)
                        .name(fusion.GeneNameDown)
                        .transcript(fusion.TranscriptUp)
                        .namePosition(0)
                        .build() :
                ImmutableGene.builder()
                        .type(finalDownGene.AnnotationType)
                        .chromosome(finalDownGene.Chromosome)
                        .start(fusion.PosDown)
                        .end(finalDownGene.ExonEnd)
                        .strand(fusion.StrandDown)
                        .name(fusion.GeneNameDown)
                        .transcript(fusion.TranscriptDown)
                        .namePosition(0)
                        .build();
    }

    static GenomeRegion convertRegion(int strand, final GenomeRegion reference, final GenomeRegion region)
    {
        final int start;
        final int end;
        if(strand < 0)
        {
            start = reference.end() - region.end();
            end = reference.end() - region.start();
        }
        else
        {
            start = region.start() - reference.start();
            end = Math.min(reference.end(), region.end()) - reference.start();
        }

        return GenomeRegions.create(region.chromosome(), Math.max(0, start), end);
    }
}
