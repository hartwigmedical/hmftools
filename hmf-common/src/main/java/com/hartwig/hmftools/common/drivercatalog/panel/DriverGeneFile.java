package com.hartwig.hmftools.common.drivercatalog.panel;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class DriverGeneFile {

    private static final String DELIMITER = "\t";

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").
                add("gene")
                .add("deletionBand")
                .add("reportMissense")
                .add("reportTruncation")
                .add("reportSplice")
                .add("reportDisruption")
                .add("reportAmplification")
                .add("favorMultiHitAndBiallelic")
                .toString();
    }

    @NotNull
    private static String toString(final DriverGene gene) {
        return new StringJoiner(DELIMITER).add(gene.gene())
                .add("\"" + (gene.deletionBand() == null ? "" : gene.deletionBand()) + "\"")
                .add(String.valueOf(gene.reportMissense()))
                .add(String.valueOf(gene.reportTruncation()))
                .add(String.valueOf(gene.reportSplice()))
                .add(String.valueOf(gene.reportDisruption()))
                .add(String.valueOf(gene.reportAmplification()))
                .add(String.valueOf(gene.favorMultiHitAndBiallelic()))
                .toString();
    }

    @NotNull
    private static DriverGene fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        ImmutableDriverGene.Builder builder = ImmutableDriverGene.builder()
                .gene(values[0])
                .deletionBand(values[1])
                .reportMissense(Boolean.parseBoolean(values[2]))
                .reportTruncation(Boolean.parseBoolean(values[3]))
                .reportSplice(Boolean.parseBoolean(values[4]))
                .reportDisruption(Boolean.parseBoolean(values[5]))
                .reportAmplification(Boolean.parseBoolean(values[6]))
                .favorMultiHitAndBiallelic(Boolean.parseBoolean(values[7]));
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
