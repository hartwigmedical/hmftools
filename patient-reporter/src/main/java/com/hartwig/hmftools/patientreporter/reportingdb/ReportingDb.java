package com.hartwig.hmftools.patientreporter.reportingdb;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.lims.LimsCoreCohort;
import com.hartwig.hmftools.common.lims.LimsStudy;
import com.hartwig.hmftools.common.utils.io.reader.LineReader;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ReportingDb {

    private static final Logger LOGGER = LogManager.getLogger(ReportingDb.class);

    private static final String DELIMITER = "\t";
    private static final String NA_STRING = "N/A";

    private ReportingDb() {
    }

    public static void addSequenceReportToReportingDb(@NotNull String reportingDbTsv, @NotNull AnalysedPatientReport report)
            throws IOException {
        String sampleId = report.sampleReport().tumorSampleId();

        LimsStudy study = LimsStudy.fromSampleId(sampleId);

        if (requiresSummary(sampleId, study) && report.clinicalSummary().isEmpty()) {
            LOGGER.warn("Skipping addition to reporting db, missing summary for sample '{}'!", sampleId);
        } else if (study != LimsStudy.NON_CANCER_STUDY) {
            String tumorBarcode = report.sampleReport().tumorSampleBarcode();
            String cohort = !report.sampleReport().cohort().isEmpty() ? report.sampleReport().cohort() : NA_STRING;
            String reportDate = ReportResources.REPORT_DATE;
            String purity = new DecimalFormat("0.00").format(report.impliedPurity());

            boolean hasReliableQuality = report.hasReliableQuality();
            boolean hasReliablePurity = report.hasReliablePurity();

            String reportType = Strings.EMPTY;
            if (report.isCorrectedReport()) {
                if (hasReliablePurity && Integer.valueOf(purity) >= 0.20) {
                    reportType = "dna_analysis_report_corrected";
                } else {
                    reportType = "dna_analysis_report_low_corrected";
                }
            } else {
                if (hasReliablePurity && Integer.valueOf(purity) >= 0.20) {
                    reportType = "dna_analysis_report";
                } else {
                    reportType = "dna_analysis_report_low";
                }
            }

            addToReportingDb(reportingDbTsv,
                    tumorBarcode,
                    sampleId,
                    cohort,
                    reportType,
                    reportDate,
                    purity,
                    hasReliableQuality,
                    hasReliablePurity);
        }
    }

    @VisibleForTesting
    static boolean requiresSummary(@NotNull String sampleId, @NotNull LimsStudy study) {
        LimsCoreCohort coreCohort = LimsCoreCohort.fromSampleId(sampleId);

        if (study == LimsStudy.WIDE) {
            return true;
        } else {
            return study == LimsStudy.CORE && coreCohort != LimsCoreCohort.CORELR02 && coreCohort != LimsCoreCohort.CORERI02;
        }
    }

    private static void addToReportingDb(@NotNull String reportingDbTsv, @NotNull String tumorBarcode, @NotNull String sampleId,
            @NotNull String cohort, @NotNull String reportType, @NotNull String reportDate, @NotNull String purity,
            boolean hasReliableQuality, boolean hasReliablePurity) throws IOException {
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
                    tumorBarcode + "\t" + sampleId + "\t" + cohort + "\t" + reportDate + "\t" + reportType + "\t" + purity + "\t"
                            + hasReliableQuality + "\t" + hasReliablePurity + "\n";
            appendToTsv(reportingDbTsv, stringToAppend);
        }
    }

    public static void addQCFailReportToReportingDb(@NotNull String reportingDbTsv, @NotNull QCFailReport report) throws IOException {
        String sampleId = report.sampleReport().tumorSampleId();
        String cohort = !report.sampleReport().cohort().isEmpty() ? report.sampleReport().cohort() : NA_STRING;
        String tumorBarcode = report.sampleReport().tumorSampleBarcode();
        String reportDate = ReportResources.REPORT_DATE;

        String reportType = report.isCorrectedReport() ? report.reason().identifier() + "_corrected" : report.reason().identifier();

        LimsStudy study = LimsStudy.fromSampleId(sampleId);
        if (study != LimsStudy.NON_CANCER_STUDY) {
            boolean present = false;
            for (ReportingEntry entry : read(reportingDbTsv)) {
                if (!present && sampleId.equals(entry.sampleId()) && tumorBarcode.equals(entry.tumorBarcode())
                        && reportType.equals(entry.reportType()) && reportDate.equals(entry.reportDate())) {
                    LOGGER.warn("Sample {} has already been reported with report type '{}' on {}!", sampleId, reportType, reportDate);
                    present = true;
                }
            }

            if (!present) {
                LOGGER.info("Adding {} to reporting db at {} with type '{}'", sampleId, reportingDbTsv, reportType);
                String stringToAppend =
                        tumorBarcode + "\t" + sampleId + "\t" + cohort + "\t" + reportDate + "\t" + reportType + "\t" + NA_STRING + "\t"
                                + NA_STRING + "\t" + NA_STRING + "\n";
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

    private static void appendToTsv(@NotNull String reportDatesTsv, @NotNull String stringToAppend) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(reportDatesTsv, true));
        writer.write(stringToAppend);
        writer.close();
    }
}
