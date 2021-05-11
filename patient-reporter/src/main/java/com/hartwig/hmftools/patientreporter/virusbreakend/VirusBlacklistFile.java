package com.hartwig.hmftools.patientreporter.virusbreakend;

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

public class VirusBlacklistFile {

    private static final Logger LOGGER = LogManager.getLogger(VirusDbFile.class);

    private static final String SEPARATOR = "\t";

    private VirusBlacklistFile() {
    }


    @NotNull
    public static VirusBlackListModel buildFromTsv(@NotNull String virusBlacklistTsv) throws IOException {
        List<String> linesVirusBlacklist = Files.readAllLines(new File(virusBlacklistTsv).toPath());

        Map<Integer, String> virusBlacklistToMap = Maps.newHashMap();

        for (String line : linesVirusBlacklist.subList(1, linesVirusBlacklist.size())) {
            String[] parts = line.split(SEPARATOR);
            if (parts.length == 3) {
                int id = Integer.parseInt(parts[0].trim());
                String taxidtype = parts[1].trim();
                virusBlacklistToMap.put(id, taxidtype);
            } else {
                LOGGER.warn("Suspicious line detected in virus db tsv: {}", line);
            }
        }

        return new VirusBlackListModel(virusBlacklistToMap);
    }
}
