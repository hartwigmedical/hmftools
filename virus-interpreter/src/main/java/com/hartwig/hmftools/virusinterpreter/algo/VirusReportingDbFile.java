package com.hartwig.hmftools.virusinterpreter.algo;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class VirusReportingDbFile {

    private static final Logger LOGGER = LogManager.getLogger(VirusReportingDbFile.class);

    private static final String SEPARATOR = "\t";

    private VirusReportingDbFile() {
    }

    @NotNull
    public static VirusReportingDbModel buildFromTsv(@NotNull String virusReportingDbTsv) throws IOException {
        List<String> linesVirusReportingDb = Files.readAllLines(new File(virusReportingDbTsv).toPath());

        Map<Integer, VirusReportingDb> speciesVirusReportingDbMap = Maps.newHashMap();

        for (String line : linesVirusReportingDb.subList(1, linesVirusReportingDb.size())) {
            String[] parts = line.split(SEPARATOR);
            ImmutableVirusReportingDb.Builder virusReportingDbBuilder = ImmutableVirusReportingDb.builder();
            if (parts.length == 6) {
                int speciesTaxid = Integer.parseInt(parts[0].trim());
                virusReportingDbBuilder.virusInterpretation(parts[1].trim())
                        .integratedMinimalCoverage(parts[2].trim().isEmpty() ? null : Integer.parseInt(parts[2].trim()))
                        .nonIntegratedMinimalCoverage(parts[3].trim().isEmpty() ? null : Integer.parseInt(parts[3].trim()))
                        .virusDriverLikelihoodType(VirusLikelihoodType.valueOf(parts[4]));
                speciesVirusReportingDbMap.put(speciesTaxid, virusReportingDbBuilder.build());
            } else {
                LOGGER.warn("Suspicious line detected in virus reporting db tsv: {}", line);
            }
        }

        return new VirusReportingDbModel(speciesVirusReportingDbMap);
    }
}
