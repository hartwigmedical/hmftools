package com.hartwig.hmftools.patientreporter.virusbreakend;

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

public final class VirusInterpretationFile {

    private static final Logger LOGGER = LogManager.getLogger(VirusInterpretationFile.class);

    private static final String SEPARATOR = "\t";

    private VirusInterpretationFile() {
    }

    @NotNull
    public static VirusInterpretationModel buildFromTsv(@NotNull String virusInterpretationTsv) throws IOException {
        List<String> linesVirusInterpretation = Files.readAllLines(new File(virusInterpretationTsv).toPath());

        Map<Integer, VirusInterpretation> speciesToInterpretationMap = Maps.newHashMap();

        for (String line : linesVirusInterpretation.subList(1, linesVirusInterpretation.size())) {
            String[] parts = line.split(SEPARATOR);
            if (parts.length == 2) {
                int speciesTaxid = Integer.parseInt(parts[0].trim());
                VirusInterpretation interpretation = VirusInterpretation.valueOf(parts[1].trim());
                speciesToInterpretationMap.put(speciesTaxid, interpretation);
            } else {
                LOGGER.warn("Suspicious line detected in virus interpretation tsv: {}", line);
            }
        }

        return new VirusInterpretationModel(speciesToInterpretationMap);
    }
}
