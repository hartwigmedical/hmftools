package com.hartwig.hmftools.virusinterpreter.algo;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class VirusWhitelistFile {

    private static final Logger LOGGER = LogManager.getLogger(VirusWhitelistFile.class);

    private static final String SEPARATOR = "\t";

    private VirusWhitelistFile() {
    }

    @NotNull
    public static VirusWhitelistModel buildFromTsv(@NotNull String virusInterpretationTsv) throws IOException {
        List<String> linesVirusWhiteList = Files.readAllLines(new File(virusInterpretationTsv).toPath());

        Map<Integer, VirusWhitelist> speciesVirusWhitelistMap = Maps.newHashMap();

        for (String line : linesVirusWhiteList.subList(1, linesVirusWhiteList.size())) {
            String[] parts = line.split(SEPARATOR);
            if (parts.length == 5) {
                int speciesTaxid = Integer.parseInt(parts[0].trim());
                VirusWhitelist virusWhitelist = ImmutableVirusWhitelist.builder()
                        .taxidSpecies(speciesTaxid)
                        .reportOnSummary(Boolean.parseBoolean(parts[1].trim()))
                        .virusInterpretation(parts[2].trim())
                        .integratedMinimalCoverage(parts[3].trim().equals(Strings.EMPTY) ? null : Integer.parseInt(parts[3].trim()))
                        .nonintegratedMinimalCoverage(parts[4].trim().equals(Strings.EMPTY) ? null : Integer.parseInt(parts[4].trim()))
                        .build();
                speciesVirusWhitelistMap.put(speciesTaxid, virusWhitelist);
            } else {
                LOGGER.warn("Suspicious line detected in virus interpretation tsv: {}", line);
            }
        }

        return new VirusWhitelistModel(speciesVirusWhitelistMap);
    }
}
