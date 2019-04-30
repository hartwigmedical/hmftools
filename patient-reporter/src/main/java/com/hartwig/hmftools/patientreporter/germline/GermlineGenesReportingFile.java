package com.hartwig.hmftools.patientreporter.germline;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.io.reader.LineReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class GermlineGenesReportingFile {
    private static final Logger LOGGER = LogManager.getLogger(GermlineGenesReportingFile.class);
    private static final String SEPARATOR = ",";

    private GermlineGenesReportingFile() {
    }

    @NotNull
    public static GermlineGenesReporting buildFromCsv(@NotNull String germlineGenesCsv) throws IOException {
        List<String> linesGermlineGenes = LineReader.build().readLines(new File(germlineGenesCsv).toPath(), line -> line.length() > 0);

        Map<String, Boolean> germlineGenesMap = Maps.newHashMap();

        for (String line : linesGermlineGenes) {
            String[] parts = line.split(SEPARATOR);

            if (parts.length == 2) {
                String gene = parts[0].trim();
                String classificationGene = parts[1].trim().toLowerCase();
                switch (classificationGene) {
                    case "true":
                        germlineGenesMap.put(gene, true);
                        break;
                    case "false":
                        germlineGenesMap.put(gene, false);
                        break;
                    default:
                        LOGGER.warn("Could not interpret classification in germline reporting genes: " + classificationGene);
                }
            } else {
                LOGGER.warn("Suspicious line detected in germline reporting genes: " + line);
            }
        }

        return ImmutableGermlineGenesReporting.of(germlineGenesMap);
    }
}


