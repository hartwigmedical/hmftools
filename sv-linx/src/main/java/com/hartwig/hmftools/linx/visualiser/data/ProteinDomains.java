package com.hartwig.hmftools.linx.visualiser.data;

import static com.hartwig.hmftools.linx.visualiser.data.FusedExons.convertRegion;

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
    public static final String UTR = "UTR/Non-coding";

    public static List<ProteinDomain> proteinDomainsInGenes(@NotNull final List<Gene> genes,
            @NotNull final List<ProteinDomain> proteinDomains)
    {
        final Predicate<ProteinDomain> matchesGene = domain -> genes.stream().anyMatch(gene -> gene.overlaps(domain));
        return proteinDomains.stream().filter(matchesGene).collect(Collectors.toList());
    }

    @NotNull
    public static List<ProteinDomain> proteinDomainsInFusion(@NotNull final Fusion fusion, @NotNull final List<FusedExon> fusedExons,
            @NotNull final List<ProteinDomain> proteinDomains)
    {
        final List<ProteinDomain> result = Lists.newArrayList();
        if (fusedExons.isEmpty())
        {
            return result;
        }

        final FusedExon firstUpExon = fusedExons.get(0);
        final GenomeRegion upGeneRegion = upGeneRegion(fusion, firstUpExon);

        final FusedExon finalDownExon = fusedExons.get(fusedExons.size() - 1);
        final GenomeRegion downGeneRegion = downGeneRegion(fusion, finalDownExon);

        for (ProteinDomain unadjustedDomain : proteinDomains)
        {
            if (unadjustedDomain.overlaps(upGeneRegion))
            {
                final GenomeRegion convertedDomain = convertRegion(fusion.strandUp(), upGeneRegion, unadjustedDomain);
                final ProteinDomain domain = ImmutableProteinDomain.builder().from(unadjustedDomain)
                        .chromosome(fusion.name())
                        .start(Math.max(convertedDomain.start(), firstUpExon.geneStart()))
                        .end(Math.min(convertedDomain.end(), firstUpExon.geneEnd()))
                        .build();

                result.add(domain);
            }

            if (unadjustedDomain.overlaps(downGeneRegion))
            {
                final GenomeRegion convertedDomain = convertRegion(fusion.strandDown(), downGeneRegion, unadjustedDomain);
                final ProteinDomain domain = ImmutableProteinDomain.builder().from(unadjustedDomain)
                        .chromosome(fusion.name())
                        .start(Math.max(convertedDomain.start() + firstUpExon.geneEnd(), finalDownExon.geneStart()))
                        .end(Math.min(convertedDomain.end() + firstUpExon.geneEnd(), finalDownExon.geneEnd()))
                        .build();

                result.add(domain);
            }
        }

        return result;
    }

    @NotNull
    private static GenomeRegion upGeneRegion(@NotNull final Fusion fusion, @NotNull final FusedExon firstUpGene)
    {
        final long upGeneLength = firstUpGene.geneEnd() - firstUpGene.geneStart();
        final long upGeneStart = fusion.strandUp() < 0 ? fusion.positionUp() : fusion.positionUp() - upGeneLength;

        return GenomeRegions.create(fusion.chromosomeUp(), upGeneStart, upGeneStart + upGeneLength);
    }

    @NotNull
    private static GenomeRegion downGeneRegion(@NotNull final Fusion fusion, @NotNull final FusedExon finalDownExon)
    {
        final long downGeneLength = finalDownExon.geneEnd() - finalDownExon.geneStart();
        final long downGeneStart = fusion.strandDown() < 0 ? fusion.positionDown() - downGeneLength : fusion.positionDown();

        return GenomeRegions.create(fusion.chromosomeDown(), downGeneStart, downGeneStart + downGeneLength);
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
                .name(utr(file.Info))
                .transcript(file.Transcript)
                .build();
    }

    @NotNull
    private static String utr(@NotNull final String utr) {
        switch (utr) {
            case "Non Coding":
            case "5-Prime UTR":
            case "3-Prime UTR":
                return UTR;
        }

        return utr;
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
                .add("transcript")
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
                .add(domain.transcript())
                .toString();
    }

    private static String hexColor(@NotNull Color color)
    {
        return String.format("#%02X%02X%02X", color.getRed(), color.getGreen(), color.getBlue());
    }

}
