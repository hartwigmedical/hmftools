package com.hartwig.hmftools.common.variant.structural.linx;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.common.fusion.ImmutableReportableGeneFusion;
import com.hartwig.hmftools.common.fusion.ReportableGeneFusion;

import org.jetbrains.annotations.NotNull;

public final class ReportableGeneFusionFile {

    private static final String DELIMITER = "\t";

    @NotNull
    public static List<ReportableGeneFusion> read(final String filePath) throws IOException {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    @NotNull
    private static List<ReportableGeneFusion> fromLines(@NotNull List<String> lines) {
        return lines.stream().filter(x -> !x.startsWith("geneStart")).map(ReportableGeneFusionFile::fromString).collect(toList());
    }

    @NotNull
    private static ReportableGeneFusion fromString(@NotNull String line) {
        String[] values = line.split(DELIMITER);

        int index = 0;

        return ImmutableReportableGeneFusion.builder()
                .geneStart(values[index++])
                .geneContextStart(values[index++])
                .geneTranscriptStart(values[index++])
                .geneEnd(values[index++])
                .geneContextEnd(values[index++])
                .geneTranscriptEnd(values[index++])
                .junctionCopyNumber(Double.valueOf(values[index]))
                .build();
    }
}
