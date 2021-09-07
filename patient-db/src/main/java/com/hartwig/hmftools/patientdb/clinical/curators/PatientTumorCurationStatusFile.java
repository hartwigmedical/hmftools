package com.hartwig.hmftools.patientdb.clinical.curators;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class PatientTumorCurationStatusFile {

    private static final String DELIMITER = "\t";

    private PatientTumorCurationStatusFile() {
    }

    public static void write(@NotNull String outputPath, @NotNull Map<String, PatientTumorCurationStatus> patientTumorCurationStatusMap)
            throws IOException {
        Files.write(new File(outputPath).toPath(), toLines(patientTumorCurationStatusMap));
    }

    @NotNull
    private static List<String> toLines(@NotNull Map<String, PatientTumorCurationStatus> patientTumorCurationStatusMap) {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        for (String patient : patientTumorCurationStatusMap.keySet()) {
            lines.add(toLine(patient, patientTumorCurationStatusMap.get(patient)));
        }
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER).add("patient").add("status").toString();
    }

    @NotNull
    private static String toLine(@NotNull String patient, @NotNull PatientTumorCurationStatus status) {
        return new StringJoiner(DELIMITER).add(patient).add(status.toString()).toString();
    }
}
