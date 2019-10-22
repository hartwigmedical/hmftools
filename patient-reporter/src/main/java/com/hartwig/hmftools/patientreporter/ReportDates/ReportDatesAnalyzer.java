package com.hartwig.hmftools.patientreporter.ReportDates;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.io.reader.LineReader;

import org.jetbrains.annotations.NotNull;

public final class ReportDatesAnalyzer {
    private static final String DELIMITER = "\t";

    private ReportDatesAnalyzer() {

    }

    @NotNull
    public static List<ReportDates> read(@NotNull String filePath) throws IOException {

        List<String> linesReportDates = LineReader.build().readLines(new File(filePath).toPath(), line -> line.length() > 0);
        List<ReportDates> reportDatesList = Lists.newArrayList();
        for (String line : linesReportDates) {
            String[] values = line.split(DELIMITER);

            if (values.length == 4) {
                reportDatesList.add(ImmutableReportDates.builder()
                        .sampleId(values[0])
                        .tumorBarcode(values[1])
                        .reportDate(values[2])
                        .sourceReport(values[3])
                        .build());
            } else {
                reportDatesList.add(ImmutableReportDates.builder()
                        .sampleId(values[0])
                        .tumorBarcode(values[1])
                        .reportDate(values[2])
                        .sourceReport(values[3])
                        .purity(values[4])
                        .status(values[5])
                        .qcStatus(values[6])
                        .build());
            }
        }

        return reportDatesList;

    }
}
