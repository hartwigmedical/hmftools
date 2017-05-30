package com.hartwig.hmftools.patientreporter.genePanel;

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

public class GenePanelModel {
    private static final Logger LOGGER = LogManager.getLogger(GenePanelModel.class);

    @NotNull
    private final Map<String, GenePanelData> dataPerGene = Maps.newHashMap();

    public GenePanelModel(@NotNull final String genePanelCsv) throws IOException, HartwigException {
        final List<String> lines = LineReader.build().readLines(new File(genePanelCsv).toPath(),
                line -> line.length() > 0);
        lines.forEach(line -> {
            final String[] columns = line.split(",");
            final String gene = columns[0];
            final String entrezId = columns[1];
            final String chromosomeBand = columns[2];
            final String role = columns.length > 3 ? columns[3] : "";
            final GenePanelData genePanelData = new GenePanelData(entrezId, role, chromosomeBand);
            addGene(gene, genePanelData);
        });
    }

    @NotNull
    public String entrezId(@NotNull final String gene) {
        return dataPerGene.get(gene).entrezId();
    }

    @NotNull
    public String chromosomeBand(@NotNull final String gene) {
        return dataPerGene.get(gene).chromosomeBand();
    }

    @NotNull
    public String type(@NotNull final String gene) {
        return dataPerGene.get(gene).type();
    }

    private void addGene(@NotNull final String gene, @NotNull final GenePanelData genePanelData) {
        if (dataPerGene.containsKey(gene)) {
            LOGGER.warn("Gene " + gene + " already present in cosmic map.");
        }
        dataPerGene.put(gene, genePanelData);
    }
}
