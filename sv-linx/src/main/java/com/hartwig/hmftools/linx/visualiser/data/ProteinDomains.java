package com.hartwig.hmftools.linx.visualiser.data;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegions;
import com.hartwig.hmftools.linx.visualiser.circos.ProteinDomainColors;
import com.hartwig.hmftools.linx.visualiser.file.VisProteinDomainFile;

import org.jetbrains.annotations.NotNull;

public class ProteinDomains
{

    private static final String DELIMITER = "\t";

    public static List<ProteinDomain> proteinDomainsInGenes(@NotNull final List<Gene> genes, @NotNull final List<ProteinDomain> proteinDomains)
    {
        final Predicate<ProteinDomain> matchesGene = domain -> genes.stream().anyMatch(gene -> gene.overlaps(domain));
        return proteinDomains.stream().filter(matchesGene).collect(Collectors.toList());
    }

    @NotNull
    public static List<ProteinDomain> proteinDomainsInFusion(@NotNull final Fusion fusion, @NotNull final List<FusedExon> fusedExons, @NotNull final List<ProteinDomain> proteinDomains)
    {
        final List<ProteinDomain> result = Lists.newArrayList();
        if (fusedExons.isEmpty())
        {
            return result;
        }

        final FusedExon firstUpGene = fusedExons.get(0);
        final long upGeneStart = firstUpGene.unadjustedGeneStart();
        final long upGeneEnd = fusion.positionUp();

        final FusedExon firstDownExon = fusedExons.stream().filter(x -> x.gene().equals(fusion.geneDown())).findFirst().get();
        final FusedExon finalDownExon = fusedExons.get(fusedExons.size() - 1);
        final long downGeneStart = fusion.positionDown();
        final long downGeneEnd = downGeneStart + finalDownExon.end() - firstDownExon.start();

        final GenomeRegion upGeneRegion =
                GenomeRegions.create(fusion.chromosomeUp(), Math.min(upGeneStart, upGeneEnd), Math.max(upGeneStart, upGeneEnd));
        final GenomeRegion downGeneRegion =
                GenomeRegions.create(fusion.chromosomeDown(), Math.min(downGeneStart, downGeneEnd), Math.max(downGeneStart, downGeneEnd));

        final long additionalDownOffset = finalDownExon.geneStart();

        for (ProteinDomain unadjustedDomain : proteinDomains)
        {
            if (unadjustedDomain.overlaps(upGeneRegion))
            {
                long unconstrainedStart = start(fusion.strandUp(), upGeneStart, unadjustedDomain);
                long unconstrainedEnd = end(fusion.strandUp(), upGeneStart, unadjustedDomain);

                ProteinDomain domain = ImmutableProteinDomain.builder().from(unadjustedDomain)
                        .chromosome(fusion.name())
                        .start(Math.max(unconstrainedStart, firstUpGene.geneStart()))
                        .end(Math.min(unconstrainedEnd, firstUpGene.geneEnd()))
                        .build();

                result.add(domain);
            }

            if (unadjustedDomain.overlaps(downGeneRegion))
            {
                long unconstrainedStart = start(fusion.strandDown(), downGeneStart, unadjustedDomain) + additionalDownOffset;
                long unconstrainedEnd = end(fusion.strandDown(), downGeneStart, unadjustedDomain) + additionalDownOffset;

                ProteinDomain domain = ImmutableProteinDomain.builder().from(unadjustedDomain)
                        .chromosome(fusion.name())
                        .start(Math.max(unconstrainedStart, finalDownExon.geneStart()))
                        .end(Math.min(unconstrainedEnd, finalDownExon.geneEnd()))
                        .build();

                result.add(domain);
            }

        }

        return result;
    }

    @NotNull
    public static List<ProteinDomain> readProteinDomains(@NotNull final String fileName) throws IOException
    {
        return VisProteinDomainFile.read(fileName).stream().map(ProteinDomains::fromFile).collect(Collectors.toList());
    }

    @NotNull
    private static ProteinDomain fromFile(@NotNull final VisProteinDomainFile file)
    {
        return ImmutableProteinDomain.builder()
                .sampleId(file.SampleId)
                .clusterId(file.ClusterId)
                .chromosome(file.Chromosome)
                .start(file.Start)
                .end(file.End)
                .name(file.Info)
                .build();

    }

    private static long start(int strand, long offset, GenomeRegion region)
    {
        return strand < 0 ? offset - region.end() : region.start() - offset;
    }

    private static long end(int strand, long offset, GenomeRegion region)
    {
        return strand < 0 ? offset - region.start() : region.end() - offset;
    }

    public static void write(@NotNull final String fileName, @NotNull final ProteinDomainColors colors,
            @NotNull final List<ProteinDomain> domains) throws IOException
    {
        Files.write(new File(fileName).toPath(), toLines(colors, domains));
    }

    @NotNull
    static List<String> toLines(@NotNull final ProteinDomainColors colors, @NotNull final List<ProteinDomain> domains)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        domains.stream().map(x -> toString(colors, x)).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER).add("sampleId")
                .add("clusterId")
                .add("chromosome")
                .add("start")
                .add("end")
                .add("name")
                .add("color")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final ProteinDomainColors colors, @NotNull final ProteinDomain domain)
    {
        return new StringJoiner(DELIMITER)
                .add(domain.sampleId())
                .add(String.valueOf(domain.clusterId()))
                .add(String.valueOf(domain.chromosome()))
                .add(String.valueOf(domain.start()))
                .add(String.valueOf(domain.end()))
                .add(String.valueOf(domain.name()))
                .add(hexColor(colors.color(domain.name())))
                .toString();
    }

    private static String hexColor(@NotNull Color color)
    {
        return String.format("#%02X%02X%02X", color.getRed(), color.getGreen(), color.getBlue());
    }

}
