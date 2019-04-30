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

        Set<String> genes = Sets.newHashSet();
        Map<String, Boolean> germlineGenesMap = Maps.newHashMap();

        for (String line : linesGermlineGenes) {
            String[] parts = line.split(SEPARATOR);

            String gene = parts[0].trim();
            genes.add(gene);

            if (parts.length == 2) {
                if (parts[1].trim().equals("true")) {
                    germlineGenesMap.put(gene, true);
                } else if (parts[1].trim().equals("false")) {
                    germlineGenesMap.put(gene, false);
                } else {
                    LOGGER.warn("No information about notify germline of gene");
                }
            }
        }

        return ImmutableGermlineGenesReporting.of(germlineGenesMap);
    }
}


