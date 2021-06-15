package com.hartwig.hmftools.common.drivercatalog.panel;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;

import org.jetbrains.annotations.NotNull;

public final class DriverGeneFile
{
    private static final String DELIMITER = "\t";

    private DriverGeneFile()
    {
    }

    public static void write(@NotNull final String filename, @NotNull final List<DriverGene> driverGenes) throws IOException
    {
        List<DriverGene> sorted = Lists.newArrayList(driverGenes);
        Files.write(new File(filename).toPath(), toLines(sorted));
    }

    @NotNull
    public static List<DriverGene> read(@NotNull final String filename) throws IOException
    {
        return Files.readAllLines(new File(filename).toPath())
                .stream()
                .filter(x -> !x.startsWith("gene") && !x.startsWith("HG") && !x.toLowerCase().startsWith("v37"))
                .map(DriverGeneFile::fromString)
                .collect(Collectors.toList());
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER).add("gene")
                .add("reportMissense")
                .add("reportNonsense")
                .add("reportSplice")
                .add("reportDeletion")
                .add("reportDisruption")
                .add("reportAmplification")
                .add("reportSomaticHotspot")
                .add("likelihoodType")
                .add("reportGermlineVariant")
                .add("reportGermlineHotspot")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final DriverGene gene)
    {
        return new StringJoiner(DELIMITER).add(gene.gene())
                .add(String.valueOf(gene.reportMissenseAndInframe()))
                .add(String.valueOf(gene.reportNonsenseAndFrameshift()))
                .add(String.valueOf(gene.reportSplice()))
                .add(String.valueOf(gene.reportDeletion()))
                .add(String.valueOf(gene.reportDisruption()))
                .add(String.valueOf(gene.reportAmplification()))
                .add(String.valueOf(gene.reportSomaticHotspot()))
                .add(String.valueOf(gene.likelihoodType()))
                .add(String.valueOf(gene.reportGermlineVariant()))
                .add(String.valueOf(gene.reportGermlineHotspot()))
                .toString();
    }

    @NotNull
    public static DriverGene fromString(@NotNull final String line)
    {
        String[] values = line.split(DELIMITER);
        ImmutableDriverGene.Builder builder = ImmutableDriverGene.builder()
                .gene(values[0])
                .reportMissenseAndInframe(Boolean.parseBoolean(values[1].toLowerCase()))
                .reportNonsenseAndFrameshift(Boolean.parseBoolean(values[2].toLowerCase()))
                .reportSplice(Boolean.parseBoolean(values[3].toLowerCase()))
                .reportDeletion(Boolean.parseBoolean(values[4].toLowerCase()))
                .reportDisruption(Boolean.parseBoolean(values[5].toLowerCase()))
                .reportAmplification(Boolean.parseBoolean(values[6].toLowerCase()))
                .reportSomaticHotspot(Boolean.parseBoolean(values[7].toLowerCase()))
                .likelihoodType(DriverCategory.valueOf(values[8]))
                .reportGermlineVariant(DriverGeneGermlineReporting.NONE)
                .reportGermlineHotspot(DriverGeneGermlineReporting.NONE);

        if(values.length == 11)
        {
            builder.reportGermlineVariant(DriverGeneGermlineReporting.valueOf(values[9].toUpperCase()))
                    .reportGermlineHotspot(DriverGeneGermlineReporting.valueOf(values[10].toUpperCase()));
        }

        return builder.build();
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<DriverGene> driverGenes)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        driverGenes.stream().map(DriverGeneFile::toString).forEach(lines::add);
        return lines;
    }
}
