package com.hartwig.hmftools.patientreporter.reportingdb;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
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

    public static void generateOutputReportDatesSeqReport(@NotNull String reportingDbTsv, @NotNull AnalysedPatientReport report)
            throws IOException {
        String sampleId = report.sampleReport().tumorSampleId();
        String tumorBarcode = report.sampleReport().tumorSampleBarcode();

        LimsSampleType type = LimsSampleType.fromSampleId(sampleId);

        String reportDate = ReportResources.REPORT_DATE;

        String purity = new DecimalFormat("#'%'").format(report.impliedPurity() * 100);
        String hasReliablePurity = String.valueOf(report.hasReliablePurity());
        String hasReliableQuality = String.valueOf(report.hasReliableQuality());

        String reasonCorrect = report.isCorrectedReport() ? "sequence_report" + "_corrected" : "sequence_report";
        String keySample = sampleId + tumorBarcode + reportDate + reasonCorrect;
        String keySample2 = sampleId + tumorBarcode + reasonCorrect + purity + hasReliablePurity + hasReliableQuality;

        boolean present = false;
        for (ReportingEntry entry : read(reportingDbTsv)) {
            String keyFile = entry.sampleId() + entry.tumorBarcode() + entry.reportDate() + entry.sourceReport();
            String keyFile2 =
                    entry.sampleId() + entry.tumorBarcode() + entry.sourceReport() + entry.purity() + entry.status() + entry.qcStatus();

            if (keySample.equals(keyFile) || keySample2.equals(keyFile2)) {
                LOGGER.warn("Sample is already reported!");
                present = true;
            } else if (sampleId.startsWith("COLO")) {
                LOGGER.warn("It is a COLO sample. This sample will not be reported!");
                present = true;
            } else if (type.equals(LimsSampleType.WIDE) && report.clinicalSummary().isEmpty()) {
                LOGGER.warn("Add summary to report for WIDE!");
                present = true;
            } else if (type.equals(LimsSampleType.CORE)) {
                if (!sampleId.startsWith("CORE01LR") && report.clinicalSummary().isEmpty()) {
                    LOGGER.warn("Add summary to report for CORE!");
                    present = true;
                } else if (!sampleId.startsWith("CORE01RI") && report.clinicalSummary().isEmpty()) {
                    LOGGER.warn("Add summary to report for CORE!");
                    present = true;
                }
            }
        }

        if (!present) {
            LOGGER.info("Writing report date to tsv file");
            String stringForFile =
                    sampleId + "\t" + tumorBarcode + "\t" + reportDate + "\t" + reasonCorrect + "\t" + purity + "\t" + hasReliablePurity
                            + "\t" + hasReliableQuality + "\n";
            writeToTSV(stringForFile, reportingDbTsv);
        }
    }

    public static void generateOutputReportDatesQCFailReport(@NotNull String reportingDbTsv, @NotNull QCFailReport report)
            throws IOException {
        String reportDate = ReportResources.REPORT_DATE;

        String sampleId = report.sampleReport().tumorSampleId();
        String tumorBarcode = report.sampleReport().tumorSampleBarcode();

        String keySample = sampleId + tumorBarcode + reportDate + report.reason();

        boolean present = false;
        for (ReportingEntry entry : read(reportingDbTsv)) {
            String keyFile = entry.sampleId() + entry.tumorBarcode() + entry.reportDate() + entry.sourceReport();
            if (keySample.equals(keyFile)) {
                LOGGER.warn("Sample {} has already been reported!", sampleId);
                present = true;
            }
        }

        if (!present) {
            LOGGER.info("Writing entry to tsv file for {}.", sampleId);
            String stringForFile = sampleId + "\t" + tumorBarcode + "\t" + reportDate + "\t" + report.reason() + "\n";
            writeToTSV(stringForFile, reportingDbTsv);
        }
    }

    @NotNull
    @VisibleForTesting
    static List<ReportingEntry> read(@NotNull String reportingDbTsv) throws IOException {
        List<String> linesReportDates = LineReader.build().readLines(new File(reportingDbTsv).toPath(), line -> line.length() > 0);
        List<ReportingEntry> reportingEntryList = Lists.newArrayList();
        for (String line : linesReportDates.subList(1, linesReportDates.size())) {
            String[] values = line.split(DELIMITER);

            if (values.length == 4) {
                reportingEntryList.add(ImmutableReportingEntry.builder()
                        .sampleId(values[0])
                        .tumorBarcode(values[1])
                        .reportDate(values[2])
                        .sourceReport(values[3])
                        .build());
            } else {
                reportingEntryList.add(ImmutableReportingEntry.builder()
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
        return reportingEntryList;
    }

    private static void writeToTSV(@NotNull String stringForFile, @NotNull String reportDatesTSV) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(reportDatesTSV, true));
        writer.write(stringForFile);
        writer.close();
    }
}
