package com.hartwig.hmftools.patientdb.clinical.consents;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortConfig;
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortModel;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortModel;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortModelFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ConsentConfigFactory {

    private static final Logger LOGGER = LogManager.getLogger(ConsentConfigFactory.class);
    private static final String DELIMITER = "\t";

    private ConsentConfigFactory() {
    }

    @NotNull
    public static Map<String, ConsentConfig> read(@NotNull String fileName) throws IOException {
        Map<String, ConsentConfig> consentConfigMap = Maps.newHashMap();
        List<String> lines = Files.readAllLines(new File(fileName).toPath());

        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(DELIMITER);

            if (parts.length == 3) {
                ConsentConfig consentConfig = ImmutableConsentConfig.builder()
                        .pifVersion(parts[0])
                        .inHMF(parts[1].equals("Yes"))
                        .outsideEU(parts[2].equals("Yes"))
                        .pif222(null)
                        .pif222Values(null)
                        .pif221(null)
                        .pif221Values(null)
                        .pif26HMF(null)
                        .pif26HMFValues(null)
                        .pif26BUG(null)
                        .pif26BUGValues(null)
                        .build();
                consentConfigMap.put(parts[0], consentConfig);
            } else if (parts.length == 7) {
                ConsentConfig consentConfig = ImmutableConsentConfig.builder()
                        .pifVersion(parts[0])
                        .inHMF(null)
                        .outsideEU(null)
                        .pif222(parts[3])
                        .pif222Values(parts[4].isEmpty() ? null : Lists.newArrayList(parts[4].split("/")))
                        .pif221(parts[5])
                        .pif221Values(parts[6].isEmpty() ? null : Lists.newArrayList(parts[6].split("/")))
                        .pif26HMF(null)
                        .pif26HMFValues(null)
                        .pif26BUG(null)
                        .pif26BUGValues(null)
                        .build();

                consentConfigMap.put(parts[0], consentConfig);
            } else if (parts.length == 11) {
                ConsentConfig consentConfig = ImmutableConsentConfig.builder()
                        .pifVersion(parts[0])
                        .inHMF(null)
                        .outsideEU(null)
                        .pif222(parts[3])
                        .pif222Values(parts[4].isEmpty() ? null : Lists.newArrayList(parts[4].split("/")))
                        .pif221(parts[5])
                        .pif221Values(parts[6].isEmpty() ? null : Lists.newArrayList(parts[6].split("/")))
                        .pif26HMF(parts[7])
                        .pif26HMFValues(parts[8].isEmpty() ? null : Lists.newArrayList(parts[8].split("/")))
                        .pif26BUG(parts[9])
                        .pif26BUGValues(parts[10].isEmpty() ? null : Lists.newArrayList(parts[10].split("/")))
                        .build();

                consentConfigMap.put(parts[0], consentConfig);
            } else {
                LOGGER.warn("Could not properly parse line in consent config tsv '{}' with length '{}'", line, parts.length);
            }
        }

        return consentConfigMap;
    }
}
