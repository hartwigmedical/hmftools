package com.hartwig.hmftools.patientreporter.cosmic;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.io.reader.LineReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CosmicCensus {
    private static final Logger LOGGER = LogManager.getLogger(CosmicCensus.class);

    @NotNull
    private final Map<String, CosmicData> dataPerGene = Maps.newHashMap();

    public CosmicCensus(@NotNull final String cosmicCsv) throws IOException, HartwigException {
        final List<String> lines = LineReader.build().readLines(new File(cosmicCsv).toPath(),
                line -> line.length() > 0);
        lines.forEach(line -> {
            final String[] columns = line.split(",");
            final String gene = columns[0];
            final String entrezId = columns[2];
            final String chromosomeBand = columns[4];
            final String role = columns[12];
            final String synonyms = columns[17];
            final CosmicData cosmicData = new CosmicData(entrezId, role, chromosomeBand);
            addGene(gene, cosmicData);
            for (final String synonym : synonyms.split(";")) {
                addGene(synonym, cosmicData);
            }
        });
    }

    @Nullable
    public CosmicData getGeneData(@NotNull final String gene) {
        return dataPerGene.get(gene);
    }

    private void addGene(@NotNull final String gene, @NotNull final CosmicData cosmicData) {
        if (dataPerGene.containsKey(gene)) {
            LOGGER.warn("Gene " + gene + " already present in cosmic map.");
        }
        dataPerGene.put(gene, cosmicData);
    }
}
