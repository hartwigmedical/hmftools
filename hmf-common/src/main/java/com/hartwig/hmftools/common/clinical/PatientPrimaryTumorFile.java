package com.hartwig.hmftools.common.clinical;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class PatientPrimaryTumorFile {

    private static final String TAB_DELIMITER = "\t";
    private static final String DOID_DELIMITER = ";";

    private PatientPrimaryTumorFile() {
    }

    @NotNull
    public static List<PatientPrimaryTumor> read(@NotNull String filePath) throws IOException {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
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
                    .primaryTumorLocation(parts[1])
                    .primaryTumorSubLocation(parts[2])
                    .primaryTumorType(parts[3])
                    .primaryTumorSubType(parts[4])
                    .primaryTumorExtraDetails(parts[5])
                    .doids(toDOIDs(parts[6]))
                    .isOverridden(Boolean.parseBoolean(parts[7]))
                    .build());
        }

        return patientPrimaryTumors;
    }

    public static void write(@NotNull String outputPath, @NotNull List<PatientPrimaryTumor> patientPrimaryTumors) throws IOException {
        Files.write(new File(outputPath).toPath(), toLines(patientPrimaryTumors));
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
                .add("primaryTumorLocation")
                .add("primaryTumorSubLocation")
                .add("primaryTumorType")
                .add("primaryTumorSubType")
                .add("primaryTumorExtraDetails")
                .add("doids")
                .add("overridden")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull PatientPrimaryTumor patientPrimaryTumor) {
        return new StringJoiner(TAB_DELIMITER).add(patientPrimaryTumor.patientIdentifier())
                .add(patientPrimaryTumor.primaryTumorLocation())
                .add(patientPrimaryTumor.primaryTumorSubLocation())
                .add(patientPrimaryTumor.primaryTumorType())
                .add(patientPrimaryTumor.primaryTumorSubType())
                .add(patientPrimaryTumor.primaryTumorExtraDetails())
                .add(fromDOIDs(patientPrimaryTumor.doids()))
                .add(String.valueOf(patientPrimaryTumor.isOverridden()))
                .toString();
    }

    @NotNull
    @VisibleForTesting
    static List<String> toDOIDs(@NotNull String doidPart) {
        return Lists.newArrayList(doidPart.split(DOID_DELIMITER));
    }

    @NotNull
    @VisibleForTesting
    static String fromDOIDs(@NotNull List<String> doids) {
        StringJoiner joiner = new StringJoiner(DOID_DELIMITER);
        for (String doid : doids) {
            joiner.add(doid);
        }
        return joiner.toString();
    }
}
