package com.hartwig.hmftools.patientreporter.reportingdb;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.io.reader.LineReader;
import com.hartwig.hmftools.common.lims.LimsSampleType;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ReportingDb {

    private static final Logger LOGGER = LogManager.getLogger(ReportingDb.class);

    private static final String DELIMITER = "\t";

    private ReportingDb() {
    }

    public static void addSequenceReportToReportingDb(@NotNull String reportingDbTsv, @NotNull AnalysedPatientReport report)
            throws IOException {
        String sampleId = report.sampleReport().tumorSampleId();
        String tumorBarcode = report.sampleReport().tumorSampleBarcode();
        String reportDate = ReportResources.REPORT_DATE;

        String purity = new DecimalFormat("#'%'").format(report.impliedPurity() * 100);
        String hasReliablePurity = String.valueOf(report.hasReliablePurity());
        String hasReliableQuality = String.valueOf(report.hasReliableQuality());

        String reportType = report.isCorrectedReport() ? "sequence_report_corrected" : "sequence_report";

        LimsSampleType type = LimsSampleType.fromSampleId(sampleId);
        if (sampleId.startsWith("COLO")) {
            LOGGER.debug("This is a COLO sample. This sample will not be included in reporting db");
        } else if (type.equals(LimsSampleType.WIDE) && report.clinicalSummary().isEmpty()) {
            LOGGER.warn("Skipping addition to reporting db, missing summary for WIDE sample {}!", sampleId);
        } else if (type.equals(LimsSampleType.CORE) && report.clinicalSummary().isEmpty() && !sampleId.startsWith("CORE01LR")
                && !sampleId.startsWith("CORE01RI")) {
            LOGGER.warn("Skipping addition to reporting db, missing summary for CORE sample {}!", sampleId);
        } else {
            boolean present = false;
            for (ReportingEntry entry : read(reportingDbTsv)) {
                if (!present && sampleId.equals(entry.sampleId()) && tumorBarcode.equals(entry.tumorBarcode())
                        && reportType.equals(entry.reportType())) {
                    LOGGER.warn("Sample has already been reported: {} with report type {}!", sampleId, reportType);
                    present = true;
                }
            }

            if (!present) {
                LOGGER.info("Adding {} to reporting db at {} with type '{}'", sampleId, reportingDbTsv, reportType);
                String stringToAppend =
                        sampleId + "\t" + tumorBarcode + "\t" + reportDate + "\t" + reportType + "\t" + purity + "\t" + hasReliablePurity
                                + "\t" + hasReliableQuality + "\n";
                appendToTsv(reportingDbTsv, stringToAppend);
            }
        }
    }

    public static void addQCFailReportToReportingDb(@NotNull String reportingDbTsv, @NotNull QCFailReport report) throws IOException {
        String sampleId = report.sampleReport().tumorSampleId();
        String tumorBarcode = report.sampleReport().tumorSampleBarcode();
        String reportDate = ReportResources.REPORT_DATE;

        String reportType = report.reason().identifier();

        if (sampleId.startsWith("COLO")) {
            LOGGER.debug("This is a COLO sample. This sample will not be included in reporting db");
        } else {
            boolean present = false;
            for (ReportingEntry entry : read(reportingDbTsv)) {
                if (!present && sampleId.equals(entry.sampleId()) && tumorBarcode.equals(entry.tumorBarcode())
                        && reportType.equals(entry.reportType()) && reportDate.equals(entry.reportDate())) {
                    LOGGER.warn("Sample has already been reported: {} with report type {}!", sampleId, reportType);
                    present = true;
                }
            }

            if (!present) {
                LOGGER.info("Adding {} to reporting db at {} with type '{}'", sampleId, reportingDbTsv, reportType);
                String stringToAppend = sampleId + "\t" + tumorBarcode + "\t" + reportDate + "\t" + reportType + "\tN/A\tN/A\tN/A\n";
                appendToTsv(reportingDbTsv, stringToAppend);
            }
        }
    }

    @NotNull
    @VisibleForTesting
    static List<ReportingEntry> read(@NotNull String reportingDbTsv) throws IOException {
        List<String> linesReportDates = LineReader.build().readLines(new File(reportingDbTsv).toPath(), line -> line.length() > 0);
        List<ReportingEntry> reportingEntryList = Lists.newArrayList();
        for (String line : linesReportDates.subList(1, linesReportDates.size())) {
            String[] values = line.split(DELIMITER);

            reportingEntryList.add(ImmutableReportingEntry.builder()
                    .sampleId(values[0])
                    .tumorBarcode(values[1])
                    .reportDate(values[2])
                    .reportType(values[3])
                    .purity(values[4])
                    .hasReliablePurity(values[5])
                    .hasReliableQuality(values[6])
                    .build());
        }
        return reportingEntryList;
    }

    private static void appendToTsv(@NotNull String reportDatesTsv, @NotNull String stringToAppend) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(reportDatesTsv, true));
        writer.write(stringToAppend);
        writer.close();
    }
}
