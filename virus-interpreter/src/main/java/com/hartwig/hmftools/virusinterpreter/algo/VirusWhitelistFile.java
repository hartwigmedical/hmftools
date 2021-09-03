package com.hartwig.hmftools.virusinterpreter.algo;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.virus.VirusInterpretation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class VirusWhitelistFile {

    private static final Logger LOGGER = LogManager.getLogger(VirusWhitelistFile.class);

    private static final String SEPARATOR = "\t";

    private VirusWhitelistFile() {
    }

    @NotNull
    public static VirusWhitelistModel buildFromTsv(@NotNull String virusInterpretationTsv) throws IOException {
        List<String> linesVirusWhiteList = Files.readAllLines(new File(virusInterpretationTsv).toPath());

        Map<Integer, VirusInterpretation> speciesToInterpretationMap = Maps.newHashMap();

        for (String line : linesVirusWhiteList.subList(1, linesVirusWhiteList.size())) {
            String[] parts = line.split(SEPARATOR);
            if (parts.length == 2) {
                int speciesTaxid = Integer.parseInt(parts[0].trim());
                VirusInterpretation interpretation = VirusInterpretation.valueOf(parts[1].trim());
                speciesToInterpretationMap.put(speciesTaxid, interpretation);
            } else {
                LOGGER.warn("Suspicious line detected in virus interpretation tsv: {}", line);
            }
        }

        return new VirusWhitelistModel(speciesToInterpretationMap);
    }
}
