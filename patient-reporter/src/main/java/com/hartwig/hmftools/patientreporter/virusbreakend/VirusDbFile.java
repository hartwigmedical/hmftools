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

public final class VirusDbFile {

    private static final Logger LOGGER = LogManager.getLogger(VirusDbFile.class);

    private static final String SEPARATOR = "\t";

    private VirusDbFile() {
    }

    @NotNull
    public static VirusDbModel buildFromTsv(@NotNull String virusDbTsv) throws IOException {
        List<String> linesVirusDb = Files.readAllLines(new File(virusDbTsv).toPath());

        Map<Integer, String> virusIdToMap = Maps.newHashMap();

        for (String line : linesVirusDb) {
            String[] parts = line.split(SEPARATOR);
            if (parts.length == 2) {
                int id = Integer.parseInt(parts[0].trim());
                String virusName = parts[1].trim();
                virusIdToMap.put(id, virusName);
            } else {
                LOGGER.warn("Suspicious line detected in virus db tsv: {}", line);
            }
        }

        return new VirusDbModel(virusIdToMap);
    }
}
