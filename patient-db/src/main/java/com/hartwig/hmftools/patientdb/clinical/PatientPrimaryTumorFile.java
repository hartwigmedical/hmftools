package com.hartwig.hmftools.patientdb.clinical;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class PatientPrimaryTumorFile {

    private static final String TAB_DELIMITER = "\t";
    private static final String STRING_DELIMITER = ";";

    private PatientPrimaryTumorFile() {
    }

    @NotNull
    public static List<PatientPrimaryTumor> read(@NotNull String filePath) throws IOException {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull String filePath, @NotNull List<PatientPrimaryTumor> patientPrimaryTumors) throws IOException {
        Files.write(new File(filePath).toPath(), toLines(patientPrimaryTumors));
    }

    @NotNull
    @VisibleForTesting
    static List<PatientPrimaryTumor> fromLines(@NotNull List<String> lines) {
        List<PatientPrimaryTumor> patientPrimaryTumors = Lists.newArrayList();
        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(TAB_DELIMITER);
            patientPrimaryTumors.add(ImmutablePatientPrimaryTumor.builder()
                    .patientIdentifier(parts[0])
                    .location(parts[1])
                    .subLocation(parts[2])
                    .type(parts[3])
                    .subType(parts[4])
                    .extraDetails(parts[5])
                    .doids(toStringList(parts[6]))
                    .snomedConceptIds(toStringList(parts[7]))
                    .isOverridden(Boolean.parseBoolean(parts[8]))
                    .build());
        }

        return patientPrimaryTumors;
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull List<PatientPrimaryTumor> patientPrimaryTumors) {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        for (PatientPrimaryTumor patientPrimaryTumor : patientPrimaryTumors) {
            lines.add(toString(patientPrimaryTumor));
        }
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(TAB_DELIMITER).add("patientIdentifier")
                .add("location")
                .add("subLocation")
                .add("type")
                .add("subType")
                .add("extraDetails")
                .add("doids")
                .add("snomedConceptIds")
                .add("overridden")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull PatientPrimaryTumor patientPrimaryTumor) {
        return new StringJoiner(TAB_DELIMITER).add(patientPrimaryTumor.patientIdentifier())
                .add(patientPrimaryTumor.location())
                .add(patientPrimaryTumor.subLocation())
                .add(patientPrimaryTumor.type())
                .add(patientPrimaryTumor.subType())
                .add(patientPrimaryTumor.extraDetails())
                .add(fromStringList(patientPrimaryTumor.doids()))
                .add(fromStringList(patientPrimaryTumor.snomedConceptIds()))
                .add(String.valueOf(patientPrimaryTumor.isOverridden()))
                .toString();
    }

    @NotNull
    @VisibleForTesting
    static List<String> toStringList(@NotNull String stringPart) {
        return Lists.newArrayList(stringPart.split(STRING_DELIMITER));
    }

    @NotNull
    @VisibleForTesting
    static String fromStringList(@NotNull List<String> strings) {
        StringJoiner joiner = new StringJoiner(STRING_DELIMITER);
        for (String string : strings) {
            joiner.add(string);
        }
        return joiner.toString();
    }
}
