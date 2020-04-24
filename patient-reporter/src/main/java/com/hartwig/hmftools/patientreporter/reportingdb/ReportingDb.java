package com.hartwig.hmftools.patientreporter.reportingdb;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.lims.LimsCohortType;
import com.hartwig.hmftools.common.lims.LimsSampleType;
import com.hartwig.hmftools.common.utils.io.reader.LineReader;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ReportingDb {

    private static final Logger LOGGER = LogManager.getLogger(ReportingDb.class);

    private static final String DELIMITER = "\t";
    private static final String NA_STRING = "N/A";

    private ReportingDb() {
    }

    public static void addSequenceReportToReportingDb(@NotNull String reportingDbTsv, @NotNull AnalysedPatientReport report)
            throws IOException {
        String tumorBarcode = report.sampleReport().tumorSampleBarcode();
        String sampleId = report.sampleReport().tumorSampleId();
        String reportDate = ReportResources.REPORT_DATE;
        String purity = new DecimalFormat("0.00").format(report.impliedPurity());

        boolean hasReliableQuality = report.hasReliableQuality();
        boolean hasReliablePurity = report.hasReliablePurity();

        String reportType = report.isCorrectedReport() ? "sequence_report_corrected" : "sequence_report";

        LimsCohortType cohortTypeCore = LimsCohortType.fromSampleId(sampleId);
        LimsSampleType type = LimsSampleType.fromSampleId(sampleId);
        if (type == LimsSampleType.WIDE && report.clinicalSummary().isEmpty()) {
            LOGGER.warn("Skipping addition to reporting db, missing summary for WIDE sample {}!", sampleId);
        } else if (type == LimsSampleType.CORE && report.clinicalSummary().isEmpty() && cohortTypeCore != LimsCohortType.CORELR02
                && cohortTypeCore != LimsCohortType.CORERI02) {
            LOGGER.warn("Skipping addition to reporting db, missing summary for CORE sample {}!", sampleId);
        } else if (type != LimsSampleType.OTHER) {
            addToReportingDb(reportingDbTsv, tumorBarcode, sampleId, reportType, reportDate, purity, hasReliableQuality, hasReliablePurity);
        }
    }

    private static void addToReportingDb(@NotNull String reportingDbTsv, @NotNull String tumorBarcode, @NotNull String sampleId,
            @NotNull String reportType, @NotNull String reportDate, @NotNull String purity, boolean hasReliableQuality,
            boolean hasReliablePurity) throws IOException {
        boolean present = false;
        for (ReportingEntry entry : read(reportingDbTsv)) {
            if (!present && sampleId.equals(entry.sampleId()) && tumorBarcode.equals(entry.tumorBarcode())
                    && reportType.equals(entry.reportType())) {
                LOGGER.warn("Sample {} has already been reported with report type '{}'!", sampleId, reportType);
                present = true;
            }
        }

        if (!present) {
            LOGGER.info("Adding {} to reporting db at {} with type '{}'", sampleId, reportingDbTsv, reportType);
            String stringToAppend =
                    tumorBarcode + "\t" + sampleId + "\t" + reportDate + "\t" + reportType + "\t" + purity + "\t" + hasReliableQuality
                            + "\t" + hasReliablePurity + "\n";
            appendToTsv(reportingDbTsv, stringToAppend);
        }
    }

    public static void addQCFailReportToReportingDb(@NotNull String reportingDbTsv, @NotNull QCFailReport report) throws IOException {
        String sampleId = report.sampleReport().tumorSampleId();
        String tumorBarcode = report.sampleReport().tumorSampleBarcode();
        String reportDate = ReportResources.REPORT_DATE;

        String reportType = report.reason().identifier();
        String reportTypeInterpret = report.isCorrectedReport() ? reportType+"_corrected" : reportType;


        LimsSampleType type = LimsSampleType.fromSampleId(sampleId);
        if (type != LimsSampleType.OTHER) {
            boolean present = false;
            for (ReportingEntry entry : read(reportingDbTsv)) {
                if (!present && sampleId.equals(entry.sampleId()) && tumorBarcode.equals(entry.tumorBarcode())
                        && reportTypeInterpret.equals(entry.reportType()) && reportDate.equals(entry.reportDate())) {
                    LOGGER.warn("Sample {} has already been reported with report type '{}' on {}!", sampleId, reportTypeInterpret, reportDate);
                    present = true;
                }
            }

            if (!present) {
                LOGGER.info("Adding {} to reporting db at {} with type '{}'", sampleId, reportingDbTsv, reportTypeInterpret);
                String stringToAppend =
                        tumorBarcode + "\t" + sampleId + "\t" + reportDate + "\t" + reportTypeInterpret + "\t" + NA_STRING + "\t" + NA_STRING + "\t"
                                + NA_STRING + "\n";
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
                    .tumorBarcode(values[0])
                    .sampleId(values[1])
                    .reportDate(values[2])
                    .reportType(values[3])
                    .purity(values[4])
                    .hasReliableQuality(values[5])
                    .hasReliablePurity(values[6])
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
