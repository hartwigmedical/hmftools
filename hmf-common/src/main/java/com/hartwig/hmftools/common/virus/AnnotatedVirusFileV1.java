package com.hartwig.hmftools.common.virus;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;

//TODO: Move to AnnotatedVirusFile when virus interpretation is validation and added to pipeline
public class AnnotatedVirusFileV1 {

    private static final String ANNOTATED_VIRUS_EXTENSION = ".virus.annotated.tsv";

    private static final String DELIMITER = "\t";

    private AnnotatedVirusFileV1() {
    }

    @NotNull
    public static String generateFileName(@NotNull String outputDir, @NotNull String sampleId) {
        return outputDir + File.separator + sampleId + ANNOTATED_VIRUS_EXTENSION;
    }

    @NotNull
    public static List<AnnotatedVirusV1> read(@NotNull String annotatedVirusTsv) throws IOException {
        return fromLines(Files.readAllLines(new File(annotatedVirusTsv).toPath()));
    }

    @VisibleForTesting
    @NotNull
    static List<AnnotatedVirusV1> fromLines(@NotNull List<String> lines) {
        return lines.stream().skip(1).map(AnnotatedVirusFileV1::fromString).collect(toList());
    }

    @NotNull
    private static AnnotatedVirusV1 fromString(@NotNull String annotatedVirus) {
        String[] values = annotatedVirus.split(DELIMITER);
        return ImmutableAnnotatedVirusV1.builder()
                .taxid(Integer.parseInt(values[0]))
                .name(values[1])
                .qcStatus(VirusBreakendQCStatus.valueOf(values[2]))
                .integrations(Integer.parseInt(values[3]))
                .interpretation(values[4])
                .reported(Boolean.parseBoolean(values[8]))
                .build();
    }
}
