package com.hartwig.hmftools.common.clinical;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class PatientMissingDoidFile {

    private static final String DELIMITER = "\t";

    private PatientMissingDoidFile() {
    }

    public static void writeMissingDoid(@NotNull String outputPath, @NotNull Map<String, String> patients) throws IOException {
        Files.write(new File(outputPath).toPath(), toLines(patients));
    }

    @NotNull
    static List<String> toLines(@NotNull Map<String, String> patients) {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        for (String patient : patients.keySet()) {
            lines.add(toLine(patient, patients.get(patient)));
        }
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER).add("patientIdentifier").add("reason").toString();
    }

    @NotNull
    private static String toLine(@NotNull String patient, @NotNull String reason) {
        return new StringJoiner(DELIMITER).add(patient).add(reason).toString();
    }
}
