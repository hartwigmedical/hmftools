package com.hartwig.hmftools.common.clinical;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class PatientTumorLocationFile {

    private static final String TAB_DELIMITER = "\t";

    private PatientTumorLocationFile() {
    }

    @NotNull
    public static List<PatientTumorLocation> readRecordsTSV(@NotNull String filePath) throws IOException {
        List<String> lines = Files.readAllLines(new File(filePath).toPath());

        List<PatientTumorLocation> patientTumorLocations = Lists.newArrayList();
        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(TAB_DELIMITER);
            String primaryTumorLocation = parts.length > 1 ? parts[1] : Strings.EMPTY;
            String cancerSubtype = parts.length > 2 ? parts[2] : Strings.EMPTY;
            patientTumorLocations.add(ImmutablePatientTumorLocation.builder()
                    .patientIdentifier(parts[0])
                    .primaryTumorLocation(primaryTumorLocation)
                    .cancerSubtype(cancerSubtype)
                    .build());
        }

        return patientTumorLocations;
    }

    public static void writeRecordsTSV(@NotNull String outputPath, @NotNull List<PatientTumorLocation> patientTumorLocations)
            throws IOException {
        Files.write(new File(outputPath).toPath(), toLines(patientTumorLocations, TAB_DELIMITER));
    }

    @NotNull
    private static List<String> toLines(@NotNull List<PatientTumorLocation> patientTumorLocations, @NotNull String delimiter) {
        List<String> lines = Lists.newArrayList();
        lines.add(header(delimiter));
        for (PatientTumorLocation patientTumorLocation : patientTumorLocations) {
            lines.add(toString(patientTumorLocation, delimiter));
        }
        return lines;
    }

    @NotNull
    private static String header(@NotNull String delimiter) {
        return new StringJoiner(delimiter, "", "").add("patientIdentifier").add("primaryTumorLocation").add("cancerSubtype").toString();
    }

    @NotNull
    private static String toString(@NotNull PatientTumorLocation patientTumorLocation, @NotNull String delimiter) {
        return new StringJoiner(delimiter).add(patientTumorLocation.patientIdentifier())
                .add(patientTumorLocation.primaryTumorLocation())
                .add(patientTumorLocation.cancerSubtype())
                .toString();
    }
}
