package com.hartwig.hmftools.common.clinical;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class PatientMissingDoidFile {
    private static final String TAB_DELIMITER = "\t";

    private PatientMissingDoidFile(){

    }

    public static void writeMissingDoid(@NotNull String outputPath, @NotNull Map<String, String> patients) throws IOException {
        Files.write(new File(outputPath).toPath(), toLinesMissingDoid(patients));
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLinesMissingDoid(@NotNull Map<String, String> patients) {
        List<String> lines = Lists.newArrayList();
        lines.add(headerMissingDoid());
        for (String patient : patients.keySet()) {
            String patientMessage = patient + TAB_DELIMITER + patients.get(patient);
            lines.add(toStringMissingDoid(patientMessage));
        }
        return lines;
    }

    @NotNull
    private static String headerMissingDoid() {
        return new StringJoiner(TAB_DELIMITER).add("patientIdentifier").add("reason")
                .toString();
    }

    @NotNull
    private static String toStringMissingDoid(@NotNull String patient) {
        return new StringJoiner(TAB_DELIMITER).add(patient)
                .toString();
    }
}
