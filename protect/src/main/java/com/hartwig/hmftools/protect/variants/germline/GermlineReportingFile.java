package com.hartwig.hmftools.protect.variants.germline;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.io.reader.LineReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class GermlineReportingFile {

    private static final Logger LOGGER = LogManager.getLogger(GermlineReportingFile.class);
    private static final String SEPARATOR = ",";

    private GermlineReportingFile() {
    }

    @NotNull
    public static GermlineReportingModel buildFromCsv(@NotNull String germlineGenesCsv) throws IOException {
        List<String> linesGermlineGenes = LineReader.build().readLines(new File(germlineGenesCsv).toPath(), line -> line.length() > 0);

        Map<String, Boolean> germlineGenesAndNotifyMap = Maps.newHashMap();

        for (String line : linesGermlineGenes) {
            String[] parts = line.split(SEPARATOR);

            if (parts.length == 2) {
                String gene = parts[0].trim();
                String notifyGene = parts[1].trim().toLowerCase();
                switch (notifyGene) {
                    case "true":
                        germlineGenesAndNotifyMap.put(gene, true);
                        break;
                    case "false":
                        germlineGenesAndNotifyMap.put(gene, false);
                        break;
                    default:
                        LOGGER.warn("Could not interpret notification string in germline reporting genes: {}", notifyGene);
                }
            } else {
                LOGGER.warn("Suspicious line detected in germline reporting genes: {}", line);
            }
        }

        return new GermlineReportingModel(germlineGenesAndNotifyMap);
    }
}


