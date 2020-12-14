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

public class LimsCohortConfigFactory {
    private static final String DELIMITER = "\t";
    private static final Logger LOGGER = LogManager.getLogger(LimsCohortConfigFactory.class);

    @NotNull
    public static LimsCohortModel read(@NotNull final String fileName) throws IOException {
        Map<String, LimsCohortConfigData> cohortConfigMap = Maps.newHashMap();
        List<String> lines = Files.readAllLines(new File(fileName).toPath());

        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(DELIMITER);
            if (parts.length == 7) {


                LimsCohortConfigData cohortConfig = ImmutableLimsCohortConfigData.builder()
                        .cohortId(parts[0])
                        .hospitalId(parts[1].equals("TRUE"))
                        .reportGermline(parts[2].equals("TRUE"))
                        .reportGermlineFlag(parts[3].equals("TRUE"))
                        .reportConclusion(parts[4].equals("TRUE"))
                        .reportViral(parts[5].equals("TRUE"))
                        .requireHospitalId(parts[6].equals("TRUE"))
                        .requireHospitalPAId(parts[7].equals("TRUE"))
                        .hospitalPersonsStudy(parts[8].equals("TRUE"))
                        .hospitalPersonsRequester(parts[9].equals("TRUE"))
                        .outputFile(parts[10].equals("TRUE"))
                        .build();

                cohortConfigMap.put(parts[0], cohortConfig);
            } else {
                LOGGER.warn("Could not properly parse line in cohort config tsv: '{}'", line);
            }
        }

        return ImmutableLimsCohortModel.builder().limsCohortMap(cohortConfigMap).build();
    }
}
