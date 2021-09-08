package com.hartwig.hmftools.virusinterpreter.algo;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class VirusReportingFile {

    private static final Logger LOGGER = LogManager.getLogger(VirusReportingFile.class);

    private static final String SEPARATOR = "\t";

    private VirusReportingFile() {
    }

    @NotNull
    public static VirusReportingModel buildFromTsv(@NotNull String virusReportingTsv) throws IOException {
        List<String> linesVirusReporting = Files.readAllLines(new File(virusReportingTsv).toPath());

        Map<Integer, VirusReporting> speciesVirusReportingMap = Maps.newHashMap();

        for (String line : linesVirusReporting.subList(1, linesVirusReporting.size())) {
            String[] parts = line.split(SEPARATOR);
            if (parts.length == 5) {
                int speciesTaxid = Integer.parseInt(parts[0].trim());
                VirusReporting virusReporting = ImmutableVirusReporting.builder()
                        .reportOnSummary(Boolean.parseBoolean(parts[1].trim()))
                        .virusInterpretation(parts[2].trim())
                        .integratedMinimalCoverage(parts[3].trim().isEmpty() ? null : Integer.parseInt(parts[3].trim()))
                        .nonIntegratedMinimalCoverage(parts[4].trim().isEmpty() ? null : Integer.parseInt(parts[4].trim()))
                        .build();
                speciesVirusReportingMap.put(speciesTaxid, virusReporting);
            } else {
                LOGGER.warn("Suspicious line detected in virus reporting tsv: {}", line);
            }
        }

        return new VirusReportingModel(speciesVirusReportingMap);
    }
}
