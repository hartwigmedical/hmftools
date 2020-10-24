package com.hartwig.hmftools.common.clinical;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PatientTumorLocationV2File {

    private static final String TAB_DELIMITER = "\t";
    private static final String DOID_DELIMITER = ";";

    private PatientTumorLocationV2File() {
    }

    @NotNull
    public static List<PatientTumorLocationV2> read(@NotNull String filePath) throws IOException {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    @NotNull
    @VisibleForTesting
    static List<PatientTumorLocationV2> fromLines(@NotNull List<String> lines) {
        List<PatientTumorLocationV2> patientTumorLocations = Lists.newArrayList();
        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(TAB_DELIMITER);
            patientTumorLocations.add(ImmutablePatientTumorLocationV2.builder()
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

    public static void write(@NotNull String outputPath, @NotNull List<PatientTumorLocationV2> patientTumorLocations) throws IOException {
        Files.write(new File(outputPath).toPath(), toLines(patientTumorLocations));
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull List<PatientTumorLocationV2> patientTumorLocations) {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        for (PatientTumorLocationV2 patientTumorLocation : patientTumorLocations) {
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
    private static String toString(@NotNull PatientTumorLocationV2 patientTumorLocation) {
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

    @Nullable
    @VisibleForTesting
    static String fromDOIDs(@Nullable List<String> doids) {
        StringJoiner joiner = new StringJoiner(DOID_DELIMITER);
        if (doids != null) {
            for (String doid : doids) {
                joiner.add(doid);
            }
        }
        return joiner.toString();
    }
}
