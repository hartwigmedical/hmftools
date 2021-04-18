package com.hartwig.hmftools.patientreporter.germline;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class GermlineReportingFile {

    private static final String SEPARATOR = "\t";

    private GermlineReportingFile() {
    }

    @NotNull
    public static GermlineReportingModel buildFromTsv(@NotNull String germlineReportingTsv) throws IOException {
        List<String> lines = Files.readAllLines(new File(germlineReportingTsv).toPath());

        List<GermlineReportingEntry> germlineReportingEntries = Lists.newArrayList();

        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(SEPARATOR);

            germlineReportingEntries.add(ImmutableGermlineReportingEntry.builder()
                    .gene(parts[0])
                    .notifyClinicalGeneticist(GermlineCondition.extractGermlineCondition(parts[1]))
                    .conditionFilter(parts.length > 3 ? parts[2] : null)
                    .build());
        }

        return new GermlineReportingModel(germlineReportingEntries);
    }
}


