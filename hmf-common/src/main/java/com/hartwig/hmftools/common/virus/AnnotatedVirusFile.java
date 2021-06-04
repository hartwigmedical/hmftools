package com.hartwig.hmftools.common.virus;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class AnnotatedVirusFile {

    private static final String DELIMITER = "\t";

    private AnnotatedVirusFile() {
    }

    @NotNull
    public static List<AnnotatedVirus> read(@NotNull String annotatedVirusTsv) throws IOException {
        return fromLines(Files.readAllLines(new File(annotatedVirusTsv).toPath()));
    }

    public static void write(@NotNull final String annotatedVirusTsv, @NotNull List<AnnotatedVirus> annotatedViruses)
            throws IOException {
        Files.write(new File(annotatedVirusTsv).toPath(), toLines(annotatedViruses));
    }

    @VisibleForTesting
    @NotNull
    static List<String> toLines(@NotNull List<AnnotatedVirus> annotatedViruses) {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        annotatedViruses.stream().map(AnnotatedVirusFile::toString).forEach(lines::add);
        return lines;
    }

    @VisibleForTesting
    @NotNull
    static List<AnnotatedVirus> fromLines(@NotNull List<String> lines) {
        return lines.stream().skip(1).map(AnnotatedVirusFile::fromString).collect(toList());
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER).add("taxid")
                .add("name")
                .add("qcStatus")
                .add("integrations")
                .add("interpretation")
                .add("reported")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull AnnotatedVirus annotatedVirus) {
        return new StringJoiner(DELIMITER).add(String.valueOf(annotatedVirus.taxid()))
                .add(annotatedVirus.name())
                .add(annotatedVirus.qcStatus().toString())
                .add(String.valueOf(annotatedVirus.integrations()))
                .add(annotatedVirus.interpretation() != null ? annotatedVirus.interpretation().toString() : Strings.EMPTY)
                .add(String.valueOf(annotatedVirus.reported()))
                .toString();
    }

    @NotNull
    private static AnnotatedVirus fromString(@NotNull String annotatedVirus) {
        String[] values = annotatedVirus.split(DELIMITER);
        return ImmutableAnnotatedVirus.builder()
                .taxid(Integer.parseInt(values[0]))
                .name(values[1])
                .qcStatus(VirusBreakendQCStatus.valueOf(values[2]))
                .integrations(Integer.parseInt(values[3]))
                .interpretation(!values[4].isEmpty() ? VirusInterpretation.valueOf(values[4]) : null)
                .reported(Boolean.parseBoolean(values[5]))
                .build();
    }
}
