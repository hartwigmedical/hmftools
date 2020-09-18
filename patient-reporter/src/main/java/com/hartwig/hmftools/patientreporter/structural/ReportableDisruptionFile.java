package com.hartwig.hmftools.patientreporter.structural;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.common.fusion.ImmutableReportableDisruption;
import com.hartwig.hmftools.common.fusion.ReportableDisruption;

import org.jetbrains.annotations.NotNull;

public final class ReportableDisruptionFile {

    private static final String DELIMITER = "\t";

    private ReportableDisruptionFile() {
    }

    @NotNull
    public static List<ReportableDisruption> read(final String filePath) throws IOException {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    @NotNull
    private static List<ReportableDisruption> fromLines(@NotNull List<String> lines) {
        return lines.stream().filter(x -> !x.startsWith("svId")).map(ReportableDisruptionFile::fromString).collect(toList());
    }

    @NotNull
    private static ReportableDisruption fromString(@NotNull String line) {
        String[] values = line.split(DELIMITER);

        int index = 0;

        return ImmutableReportableDisruption.builder()
                .svId(Integer.parseInt(values[index++]))
                .chromosome(values[index++])
                .orientation(Byte.parseByte(values[index++]))
                .strand(Integer.parseInt(values[index++]))
                .chrBand(values[index++])
                .gene(values[index++])
                .type(values[index++])
                .junctionCopyNumber(Double.valueOf(values[index++]))
                .exonUp(Integer.parseInt(values[index++]))
                .exonDown(Integer.parseInt(values[index++]))
                .undisruptedCopyNumber(Double.parseDouble(values[index]))
                .build();
    }
}
