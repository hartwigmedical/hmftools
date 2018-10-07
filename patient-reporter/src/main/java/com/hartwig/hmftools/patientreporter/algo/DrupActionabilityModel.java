package com.hartwig.hmftools.patientreporter.algo;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.io.reader.LineReader;

import org.jetbrains.annotations.NotNull;

public class DrupActionabilityModel {

    @NotNull
    private final Set<String> actionableGenes;

    public DrupActionabilityModel(@NotNull final String drupGenesCsv) throws IOException {
        final List<String> geneLines = LineReader.build().readLines(new File(drupGenesCsv).toPath(),
                line -> line.length() > 0);
        this.actionableGenes = Sets.newHashSet(geneLines);
    }

    @NotNull
    public Set<String> actionableGenes() {
        return actionableGenes;
    }
}
