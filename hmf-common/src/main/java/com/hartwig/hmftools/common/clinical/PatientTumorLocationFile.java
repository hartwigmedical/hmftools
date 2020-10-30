package com.hartwig.hmftools.common.clinical;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class PatientTumorLocationFile {

    private static final String TAB_DELIMITER = "\t";
    private static final String DOID_DELIMITER = ";";

    private PatientTumorLocationFile() {
    }

    @NotNull
    public static List<PatientTumorLocation> read(@NotNull String filePath) throws IOException {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    @NotNull
    @VisibleForTesting
    static List<PatientTumorLocation> fromLines(@NotNull List<String> lines) {
        List<PatientTumorLocation> patientTumorLocations = Lists.newArrayList();
        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(TAB_DELIMITER);
            patientTumorLocations.add(ImmutablePatientTumorLocation.builder()
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

        return patientTumorLocations;
    }

    public static void write(@NotNull String outputPath, @NotNull List<PatientTumorLocation> patientTumorLocations) throws IOException {
        Files.write(new File(outputPath).toPath(), toLines(patientTumorLocations));
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull List<PatientTumorLocation> patientTumorLocations) {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        for (PatientTumorLocation patientTumorLocation : patientTumorLocations) {
            lines.add(toString(patientTumorLocation));
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
    private static String toString(@NotNull PatientTumorLocation patientTumorLocation) {
        return new StringJoiner(TAB_DELIMITER).add(patientTumorLocation.patientIdentifier())
                .add(patientTumorLocation.primaryTumorLocation())
                .add(patientTumorLocation.primaryTumorSubLocation())
                .add(patientTumorLocation.primaryTumorType())
                .add(patientTumorLocation.primaryTumorSubType())
                .add(patientTumorLocation.primaryTumorExtraDetails())
                .add(fromDOIDs(patientTumorLocation.doids()))
                .add(String.valueOf(patientTumorLocation.isOverridden()))
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
