package com.hartwig.hmftools.patientdb.clinical.consents;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

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

            if (parts.length == 2) {
                ConsentConfig consentConfig = ImmutableConsentConfig.builder()
                        .pifVersion(parts[0])
                        .cohort(parts[1])
                        .inHMF(null)
                        .outsideEU(null)
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
            } else if (parts.length == 4) {
                ConsentConfig consentConfig = ImmutableConsentConfig.builder()
                        .pifVersion(parts[0])
                        .cohort(parts[1])
                        .inHMF(parts[2].equals("Yes"))
                        .outsideEU(parts[3].equals("Yes"))
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
            } else if (parts.length == 8) {
                ConsentConfig consentConfig = ImmutableConsentConfig.builder()
                        .pifVersion(parts[0])
                        .cohort(parts[1])
                        .inHMF(null)
                        .outsideEU(null)
                        .pif222(parts[4])
                        .pif222Values(parts[5].isEmpty() ? null : Lists.newArrayList(parts[5].split("/")))
                        .pif221(parts[6])
                        .pif221Values(parts[7].isEmpty() ? null : Lists.newArrayList(parts[7].split("/")))
                        .pif26HMF(null)
                        .pif26HMFValues(null)
                        .pif26BUG(null)
                        .pif26BUGValues(null)
                        .build();

                consentConfigMap.put(parts[0], consentConfig);
            } else if (parts.length == 12) {
                ConsentConfig consentConfig = ImmutableConsentConfig.builder()
                        .pifVersion(parts[0])
                        .cohort(parts[1])
                        .inHMF(null)
                        .outsideEU(null)
                        .pif222(parts[4])
                        .pif222Values(parts[5].isEmpty() ? null : Lists.newArrayList(parts[5].split("/")))
                        .pif221(parts[6])
                        .pif221Values(parts[7].isEmpty() ? null : Lists.newArrayList(parts[7].split("/")))
                        .pif26HMF(parts[8])
                        .pif26HMFValues(parts[9].isEmpty() ? null : Lists.newArrayList(parts[9].split("/")))
                        .pif26BUG(parts[10])
                        .pif26BUGValues(parts[11].isEmpty() ? null : Lists.newArrayList(parts[11].split("/")))
                        .build();

                consentConfigMap.put(parts[0], consentConfig);
            } else {
                LOGGER.warn("Could not properly parse line in consent config tsv '{}' with length '{}'", line, parts.length);
            }
        }

        return consentConfigMap;
    }
}
