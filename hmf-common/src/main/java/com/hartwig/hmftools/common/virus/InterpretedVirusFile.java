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

public final class InterpretedVirusFile {

    private static final String DELIMITER = "\t";

    private InterpretedVirusFile() {
    }

    @NotNull
    public static List<InterpretedVirus> read(@NotNull String interpretedVirusTsv) throws IOException {
        return fromLines(Files.readAllLines(new File(interpretedVirusTsv).toPath()));
    }

    public static void write(@NotNull final String interpretedVirusTsv, @NotNull List<InterpretedVirus> interpretedViruses)
            throws IOException {
        Files.write(new File(interpretedVirusTsv).toPath(), toLines(interpretedViruses));
    }

    @VisibleForTesting
    @NotNull
    static List<String> toLines(@NotNull List<InterpretedVirus> interpretedViruses) {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        interpretedViruses.stream().map(InterpretedVirusFile::toString).forEach(lines::add);
        return lines;
    }

    @VisibleForTesting
    @NotNull
    static List<InterpretedVirus> fromLines(@NotNull List<String> lines) {
        return lines.stream().skip(1).map(InterpretedVirusFile::fromString).collect(toList());
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
    private static String toString(@NotNull InterpretedVirus interpretedVirus) {
        return new StringJoiner(DELIMITER).add(String.valueOf(interpretedVirus.taxid()))
                .add(interpretedVirus.name())
                .add(interpretedVirus.qcStatus().toString())
                .add(String.valueOf(interpretedVirus.integrations()))
                .add(interpretedVirus.interpretation() != null ? interpretedVirus.interpretation().toString() : Strings.EMPTY)
                .add(String.valueOf(interpretedVirus.reported()))
                .toString();
    }

    @NotNull
    private static InterpretedVirus fromString(@NotNull String interpretedVirus) {
        String[] values = interpretedVirus.split(DELIMITER);
        return ImmutableInterpretedVirus.builder()
                .taxid(Integer.parseInt(values[0]))
                .name(values[1])
                .qcStatus(VirusBreakendQCStatus.valueOf(values[2]))
                .integrations(Integer.parseInt(values[3]))
                .interpretation(!values[4].isEmpty() ? VirusInterpretation.valueOf(values[4]) : null)
                .reported(Boolean.parseBoolean(values[5]))
                .build();
    }
}
