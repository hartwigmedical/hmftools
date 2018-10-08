package com.hartwig.hmftools.patientreporter.algo;

import java.io.File;
import java.io.IOException;
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

public final class DrupActionabilityModelFactory {

    private static final Logger LOGGER = LogManager.getLogger(DrupActionabilityModel.class);

    private static final String SEPARATOR = ",";
    private static final String ONCO = "onco";
    private static final String TSG = "tsg";

    private DrupActionabilityModelFactory() {
    }

    @NotNull
    public static DrupActionabilityModel buildFromCsv(@NotNull final String drupGenesCsv) throws IOException {
        List<String> lines = LineReader.build().readLines(new File(drupGenesCsv).toPath(), line -> line.length() > 0);

        Set<String> genes = Sets.newHashSet();
        Map<String, DriverCategory> geneDriverCategoryMap = Maps.newHashMap();

        for (String line : lines) {
            String[] parts = line.split(SEPARATOR);

            String gene = parts[0].trim();
            genes.add(gene);

            if (parts.length == 2) {
                String classification = parts[1].trim();
                switch (classification) {
                    case ONCO:
                        geneDriverCategoryMap.put(gene, DriverCategory.ONCO);
                        break;
                    case TSG:
                        geneDriverCategoryMap.put(gene, DriverCategory.TSG);
                        break;
                    default:
                        LOGGER.warn("Could not resolve a classification in DRUP genes csv: " + classification);
                        break;
                }
            }
        }

        return ImmutableDrupActionabilityModel.of(genes, geneDriverCategoryMap);
    }
}
