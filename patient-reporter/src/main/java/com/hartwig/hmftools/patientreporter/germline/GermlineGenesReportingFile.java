package com.hartwig.hmftools.patientreporter.germline;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.io.reader.LineReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class GermlineGenesReportingFile {
    private static final Logger LOGGER = LogManager.getLogger(GermlineGenesReportingFile.class);

    private GermlineGenesReportingFile() {
    }

    @NotNull
    public static GermlineGenesReporting buildFromCsv(@NotNull String germlineGenesCsv, @NotNull String germlineVariantsNotifyCsv)
            throws IOException {
        List<String> linesGermlineGenes = LineReader.build().readLines(new File(germlineGenesCsv).toPath(), line -> line.length() > 0);
        List<String> linesGermlineVariantNotify =
                LineReader.build().readLines(new File(germlineVariantsNotifyCsv).toPath(), line -> line.length() > 0);

        Set<String> germlineGenes = Sets.newHashSet();
        Set<String> germlineVariantNotify = Sets.newHashSet();

        for (String lineGenes : linesGermlineGenes) {
            germlineGenes.add(lineGenes);
        }

        for (String lineGenesNotify : linesGermlineVariantNotify) {
            germlineVariantNotify.add(lineGenesNotify);
        }

        LOGGER.info(germlineGenes);
        LOGGER.info(germlineVariantNotify);
        return ImmutableGermlineGenesReporting.of(germlineGenes, germlineVariantNotify);
    }
}


