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
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExon;

import org.jetbrains.annotations.NotNull;

public class FusedExons
{
    private static final String DELIMITER = "\t";

    public static void write(final String fileName, final List<FusedExon> fusedExons) throws IOException
    {
        Files.write(new File(fileName).toPath(), toLines(fusedExons));
    }

    public static List<FusedExon> fusedExons(final Fusion fusion, final List<VisGeneExon> exons)
    {
        final List<FusedExon> result = Lists.newArrayList();

        final List<VisGeneExon> fusionExons = VisExons.fusionExons(fusion, exons);
        final List<VisGeneExon> upStreamExons = sortedUpstreamExons(fusion, fusionExons);
        final List<VisGeneExon> downStreamExons = sortedDownstreamExons(fusion, fusionExons);

        if (upStreamExons.isEmpty() || downStreamExons.isEmpty())
            return result;

        final VisGeneExon firstUpExon = upStreamExons.get(0);
        final GenomeRegion upGeneRegion = upGeneRegion(fusion, firstUpExon);
        final GenomeRegion convertedUpGeneRegion = convertRegion(fusion.strandUp(), upGeneRegion, upGeneRegion);

        final ImmutableFusedExon.Builder upFusedExonBuilder = ImmutableFusedExon.builder()
                .sampleId(fusion.sampleId())
                .clusterId(fusion.clusterId())
                .fusion(fusion.name())
                .chromosome(fusion.chromosomeUp())
                .unadjustedGeneStart(upGeneRegion.start())
                .gene(fusion.geneUp())
                .geneStart(convertedUpGeneRegion.start())
                .geneEnd(convertedUpGeneRegion.end())
                .transcript(fusion.transcriptUp());

        for (final VisGeneExon exon : upStreamExons)
        {
            final GenomeRegion convertedExon = convertRegion(fusion.strandUp(), upGeneRegion, exon);

            if(upGeneRegion.overlaps(exon))
            {
                final FusedExon fusedExon = upFusedExonBuilder
                        .start(convertedExon.start())
                        .end(convertedExon.end())
                        .rank(exon.ExonRank)
                        .skipped(exon.ExonRank > fusion.fusedExonUp())
                        .build();
                result.add(fusedExon);
            }
        }

        final VisGeneExon finalDownExon = downStreamExons.get(downStreamExons.size() - 1);
        final GenomeRegion downGeneRegion = downGene(fusion, finalDownExon);
        final GenomeRegion convertedDownGeneRegion = convertRegion(fusion.strandDown(), downGeneRegion, downGeneRegion);

        final ImmutableFusedExon.Builder downFusedExonBuilder = ImmutableFusedExon.builder().from(upFusedExonBuilder.build())
                .chromosome(fusion.chromosomeDown())
                .unadjustedGeneStart(fusion.positionDown())
                .gene(fusion.geneDown())
                .geneStart(convertedDownGeneRegion.start() + convertedUpGeneRegion.end())
                .geneEnd(convertedDownGeneRegion.end() + convertedUpGeneRegion.end())
                .transcript(fusion.transcriptDown());

        boolean intronicToExonicFusion = fusion.regionTypeUp().equals("Intronic") && fusion.regionTypeDown().equals("Exonic");

        for (int i = 0; i < downStreamExons.size(); i++)
        {
            final VisGeneExon exon = downStreamExons.get(i);
            final GenomeRegion convertedExon = convertRegion(fusion.strandDown(), downGeneRegion, exon);

            if(downGeneRegion.overlaps(exon))
            {
                final FusedExon fusedExon = downFusedExonBuilder
                        .start(convertedExon.start() + convertedUpGeneRegion.end())
                        .end(convertedExon.end() + convertedUpGeneRegion.end())
                        .rank(exon.ExonRank)
                        .skipped(exon.ExonRank < fusion.fusedExonDown() || (i == 0 && intronicToExonicFusion))
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

    private static Gene upGeneRegion(final Fusion fusion, final VisGeneExon firstUpExon)
    {
        return fusion.strandUp() < 0 ?
                ImmutableGene.builder()
                        .type(firstUpExon.AnnotationType)
                        .chromosome(firstUpExon.Chromosome)
                        .start(fusion.positionUp())
                        .end(firstUpExon.ExonEnd)
                        .strand(fusion.strandUp())
                        .name(fusion.geneDown())
                        .transcript(fusion.transcriptUp())
                        .namePosition(0)
                        .build() :
                ImmutableGene.builder()
                        .type(firstUpExon.AnnotationType)
                        .chromosome(firstUpExon.Chromosome)
                        .start(firstUpExon.ExonStart)
                        .end(fusion.positionUp())
                        .strand(fusion.strandUp())
                        .name(fusion.geneDown())
                        .transcript(fusion.transcriptDown())
                        .namePosition(0)
                        .build();
    }

    static Gene downGene(final Fusion fusion, final VisGeneExon finalDownGene)
    {
        return fusion.strandDown() < 0 ?
                ImmutableGene.builder()
                        .type(finalDownGene.AnnotationType)
                        .chromosome(finalDownGene.Chromosome)
                        .start(finalDownGene.ExonStart)
                        .end(fusion.positionDown())
                        .strand(fusion.strandDown())
                        .name(fusion.geneDown())
                        .transcript(fusion.transcriptUp())
                        .namePosition(0)
                        .build() :
                ImmutableGene.builder()
                        .type(finalDownGene.AnnotationType)
                        .chromosome(finalDownGene.Chromosome)
                        .start(fusion.positionDown())
                        .end(finalDownGene.ExonEnd)
                        .strand(fusion.strandDown())
                        .name(fusion.geneDown())
                        .transcript(fusion.transcriptDown())
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
