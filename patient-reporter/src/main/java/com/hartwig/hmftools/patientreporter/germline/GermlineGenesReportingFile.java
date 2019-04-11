package com.hartwig.hmftools.patientreporter.germline;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Set;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;

public class GermlineGenesReportingFile {

    @NotNull
    private static Set<String> germlineGenes() {
        final InputStream inputStream = GermlineGenesReportingFile.class.getResourceAsStream("/germline/GermlineGenesReported.tsv");
        return new BufferedReader(new InputStreamReader(inputStream)).lines().collect(Collectors.toSet());
    }

    @NotNull
    private static Set<String> notifyGermlineGenes() {
        final InputStream inputStream = GermlineGenesReportingFile.class.getResourceAsStream("/germline/GermlineVariantNotify.tsv");
        return new BufferedReader(new InputStreamReader(inputStream)).lines().collect(Collectors.toSet());
    }

    @NotNull
    public static GermlineGenesReporting readingLines() {
        return ImmutableGermlineGenesReporting.builder().germlineGenes(germlineGenes()).germlineGenesNotify(notifyGermlineGenes()).build();
    }
}
