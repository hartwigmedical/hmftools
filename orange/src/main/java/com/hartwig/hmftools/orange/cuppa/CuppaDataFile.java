package com.hartwig.hmftools.orange.cuppa;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class CuppaDataFile {

    private static final String FIELD_DELIMITER = ",";

    private CuppaDataFile() {
    }

    @NotNull
    public static List<CuppaData> read(@NotNull String cuppaDataCsv) throws IOException {
        List<CuppaData> cuppaData = Lists.newArrayList();
        List<String> lines = Files.readAllLines(new File(cuppaDataCsv).toPath());

        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            cuppaData.add(fromLine(line));
        }
        return cuppaData;
    }

    @NotNull
    private static CuppaData fromLine(@NotNull String line) {
        String[] values = line.split(FIELD_DELIMITER);

        return ImmutableCuppaData.builder()
                .category(values[1])
                .resultType(values[2])
                .dataType(values[3])
                .value(values[4])
                .refCancerType(values[5])
                .refValue(Double.parseDouble(values[6]))
                .build();
    }
}
