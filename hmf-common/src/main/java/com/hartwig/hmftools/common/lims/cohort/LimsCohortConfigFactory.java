package com.hartwig.hmftools.common.lims.cohort;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.stream.Collectors;
import org.jetbrains.annotations.NotNull;

public class LimsCohortConfigFactory {
    private static final String DELIMITER = "\t";

    @NotNull
    public static List<LimsCohortConfigData> read(@NotNull final String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    @NotNull
    static List<LimsCohortConfigData> fromLines(@NotNull final List<String> lines) {
        return lines.stream().skip(1).map(x -> fromString(x)).collect(Collectors.toList());
    }

    @NotNull
    private static LimsCohortConfigData fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        ImmutableLimsCohortConfigData.Builder builder = ImmutableLimsCohortConfigData.builder()
                .cohortId(values[0])
                .reportGermline(values[1])
                .reportGermlineFlag(values[2])
                .reportConclusion(values[3])
                .reportViral(values[4])
                .requireHospitalId(values[5])
                .requireHospitalPAId(values[6]);

        return builder.build();
    }
}
