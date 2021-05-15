package com.hartwig.hmftools.common.lims.cohort;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class LimsCohortModelFactory {

    private static final Logger LOGGER = LogManager.getLogger(LimsCohortModelFactory.class);
    private static final String DELIMITER = "\t";

    private LimsCohortModelFactory() {
    }

    @NotNull
    public static LimsCohortModel read(@NotNull String fileName) throws IOException {
        Map<String, LimsCohortConfig> cohortConfigMap = Maps.newHashMap();
        List<String> lines = Files.readAllLines(new File(fileName).toPath());

        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(DELIMITER);
            if (parts.length == 12) {
                LimsCohortConfig cohortConfig = ImmutableLimsCohortConfig.builder()
                        .cohortId(parts[0])
                        .sampleContainsHospitalCenterId(parts[1].equals("TRUE"))
                        .reportGermline(parts[2].equals("TRUE"))
                        .reportGermlineFlag(parts[3].equals("TRUE"))
                        .reportConclusion(parts[4].equals("TRUE"))
                        .reportViral(parts[5].equals("TRUE"))
                        .reportPeach(parts[6].equals("TRUE"))
                        .requireHospitalId(parts[7].equals("TRUE"))
                        .requireHospitalPAId(parts[8].equals("TRUE"))
                        .requireHospitalPersonsStudy(parts[9].equals("TRUE"))
                        .requireHospitalPersonsRequester(parts[10].equals("TRUE"))
                        .requireAdditionalInformationForSidePanel(parts[11].equals("TRUE"))
                        .build();

                cohortConfigMap.put(parts[0], cohortConfig);
            } else {
                LOGGER.warn("Could not properly parse line in cohort config tsv: '{}'", line);
            }
        }

        return ImmutableLimsCohortModel.builder().limsCohortMap(cohortConfigMap).build();
    }
}
