package com.hartwig.hmftools.common.reportingdb;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class ReportingDatabase {

    private static final String DELIMITER = "\t";

    private ReportingDatabase() {
    }

    @NotNull
    @VisibleForTesting
    public static List<ReportingEntry> read(@NotNull String reportingDbTsv) throws IOException {
        List<String> linesReportDates = Files.readAllLines(new File(reportingDbTsv).toPath());
        List<ReportingEntry> reportingEntryList = Lists.newArrayList();

        for (String line : linesReportDates.subList(1, linesReportDates.size())) {
            String[] values = line.split(DELIMITER);

            reportingEntryList.add(ImmutableReportingEntry.builder()
                    .tumorBarcode(values[0])
                    .sampleId(values[1])
                    .cohort(values[2])
                    .reportDate(values[3])
                    .reportType(values[4])
                    .purity(values[5])
                    .hasReliableQuality(values[6])
                    .hasReliablePurity(values[7])
                    .build());
        }
        return reportingEntryList;
    }
}
