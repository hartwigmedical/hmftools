package com.hartwig.hmftools.patientreporter.virusbreakend;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class VirusSummaryfile {

    private static final Logger LOGGER = LogManager.getLogger(VirusSummaryfile.class);

    private static final String SEPARATOR = "\t";

    private VirusSummaryfile() {
    }

    @NotNull
    public static VirusSummaryModel buildFromTsv(@NotNull String virusSummaryTsv) throws IOException {
        List<String> linesVirusSummary = Files.readAllLines(new File(virusSummaryTsv).toPath());

        Map<Integer, String> virusIdToMap = Maps.newHashMap();

        for (String line : linesVirusSummary) {
            String[] parts = line.split(SEPARATOR);
            if (parts.length == 2) {
                int id = Integer.valueOf(parts[0].trim());
                String virusName = parts[1].trim();
                virusIdToMap.put(id, virusName);
            } else {
                LOGGER.warn("Suspicious line detected in virus summary tsv: {}", line);
            }
        }

        return new VirusSummaryModel(virusIdToMap);
    }
}
