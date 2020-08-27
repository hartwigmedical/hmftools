package com.hartwig.hmftools.common.drivercatalog.panel;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;

import org.jetbrains.annotations.NotNull;

public final class DriverGeneFile {

    private static final String DELIMITER = "\t";

    private DriverGeneFile() {
    }

    @NotNull
    public static List<DriverGene> read(@NotNull final String filename) throws IOException {
        return Files.readAllLines(new File(filename).toPath())
                .stream()
                .filter(x -> !x.startsWith("gene"))
                .map(DriverGeneFile::fromString)
                .collect(Collectors.toList());
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").
                add("gene")
                .add("deletionBand")
                .add("reportMissense")
                .add("reportNonsense")
                .add("reportSplice")
                .add("reportDeletion")
                .add("reportDisruption")
                .add("reportAmplification")
                .add("reportHotspot")
                .add("likelihoodType")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final DriverGene gene) {
        return new StringJoiner(DELIMITER).add(gene.gene())
                .add(gene.deletionBand())
                .add(String.valueOf(gene.reportMissenseAndInframe()))
                .add(String.valueOf(gene.reportNonsenseAndFrameshift()))
                .add(String.valueOf(gene.reportSplice()))
                .add(String.valueOf(gene.reportDeletion()))
                .add(String.valueOf(gene.reportDisruption()))
                .add(String.valueOf(gene.reportAmplification()))
                .add(String.valueOf(gene.reportHotspot()))
                .add(String.valueOf(gene.likelihoodType()))
                .toString();
    }

    @NotNull
    public static DriverGene fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        ImmutableDriverGene.Builder builder = ImmutableDriverGene.builder()
                .gene(values[0])
                .deletionBand(values[1])
                .reportMissenseAndInframe(Boolean.parseBoolean(values[2]))
                .reportNonsenseAndFrameshift(Boolean.parseBoolean(values[3]))
                .reportSplice(Boolean.parseBoolean(values[4]))
                .reportDeletion(Boolean.parseBoolean(values[5]))
                .reportDisruption(Boolean.parseBoolean(values[6]))
                .reportAmplification(Boolean.parseBoolean(values[7]))
                .reportHotspot(Boolean.parseBoolean(values[8]))
                .likelihoodType(DriverCategory.valueOf(values[9]));
        return builder.build();
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<DriverGene> driverGenes) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        driverGenes.stream().map(DriverGeneFile::toString).forEach(lines::add);
        return lines;
    }

    public static void write(@NotNull final String filename, @NotNull final List<DriverGene> driverGenes) throws IOException {
        List<DriverGene> sorted = Lists.newArrayList(driverGenes);
        Collections.sort(sorted);
        Files.write(new File(filename).toPath(), toLines(sorted));
    }
}
