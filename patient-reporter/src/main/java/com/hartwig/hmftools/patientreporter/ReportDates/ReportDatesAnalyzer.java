package com.hartwig.hmftools.patientreporter.ReportDates;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.io.reader.LineReader;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.common.purple.qc.PurpleQCFile;
import com.hartwig.hmftools.common.purple.qc.PurpleQCStatus;
import com.hartwig.hmftools.patientreporter.SampleMetadata;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ReportDatesAnalyzer {
    private static final String DELIMITER = "\t";
    private static final Logger LOGGER = LogManager.getLogger(ReportDatesAnalyzer.class);

    private ReportDatesAnalyzer() {

    }

    @NotNull
    public static List<ReportDates> read(@NotNull String filePath) throws IOException {

        List<String> linesReportDates = LineReader.build().readLines(new File(filePath).toPath(), line -> line.length() > 0);
        List<ReportDates> reportDatesList = Lists.newArrayList();
        for (String line : linesReportDates.subList(1, linesReportDates.size())) {
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

    public static void generateOutputReportDatesSeqRapports(@NotNull String reportDatesTsv, @NotNull String purplePurityTsv,
            @NotNull SampleMetadata sampleMetadata, @NotNull String purpleQCFile, boolean correctReport, String clinicalSummary,
            boolean addToFile) throws IOException {

        String reportDate = currentDate();

        String sampleId = sampleMetadata.tumorSampleId();
        String tumorBarcode = sampleMetadata.tumorSampleBarcode();

        PurityContext purityContext = FittedPurityFile.read(purplePurityTsv);
        final PurpleQC purpleQC = PurpleQCFile.read(purpleQCFile);

        String purity = new DecimalFormat("#'%'").format(purityContext.bestFit().purity() * 100);
        FittedPurityStatus status = purityContext.status();
        PurpleQCStatus qcStatus = purpleQC.status();

        List<ReportDates> allReportDates = ReportDatesAnalyzer.read(reportDatesTsv);
        String reasonCorrect = correctReport ? "sequence_report" + "_corrected" : "sequence_report";
        String keySample = sampleId + tumorBarcode + reportDate + reasonCorrect;
        String keySample2 = sampleId + tumorBarcode + reasonCorrect + purity + status + qcStatus;

        boolean present = false;
        for (ReportDates dates : allReportDates) {
            String keyFile = dates.sampleId() + dates.tumorBarcode() + dates.reportDate() + dates.sourceReport();
            String keyFile2 =
                    dates.sampleId() + dates.tumorBarcode() + dates.sourceReport() + dates.purity() + dates.status() + dates.qcStatus();

            if (keySample.equals(keyFile) || keySample2.equals(keyFile2)) {
                LOGGER.warn("Sample is already reported");
                present = true;
            } else if (sampleId.startsWith("WIDE") && clinicalSummary.isEmpty()) {
                LOGGER.warn("Add summary to report for WIDE");
                present = true;
            } else if (sampleId.startsWith("CORE")) {
                if (!sampleId.startsWith("CORE01LR") && clinicalSummary.isEmpty()) {
                    LOGGER.warn("Add summary to report for CORE");
                    present = true;
                } else if (!sampleId.startsWith("CORE01RI") && clinicalSummary.isEmpty()) {
                    LOGGER.warn("Add summary to report for CORE");
                    present = true;
                }
            }
        }

        if (!present && addToFile) {
            LOGGER.info("Writing report date to tsv file");
            String stringForFile =
                    sampleId + "\t" + tumorBarcode + "\t" + reportDate + "\t" + reasonCorrect + "\t" + purity + "\t" + status + "\t"
                            + qcStatus + "\n";
            writeToTSV(stringForFile, reportDatesTsv);
        }
    }

    public static void generateOutputReportDatesQCFailReport(QCFailReason reason, @NotNull String reportDatesTSV,
            @NotNull SampleMetadata sampleMetadata, boolean addToFile) throws IOException {
        String reportDate = currentDate();

        String sampleId = sampleMetadata.tumorSampleId();
        String tumorBarcode = sampleMetadata.tumorSampleBarcode();

        List<ReportDates> allReportDates = ReportDatesAnalyzer.read(reportDatesTSV);

        String keySample = sampleId + tumorBarcode + reportDate + reason;

        boolean present = false;
        for (ReportDates dates : allReportDates) {
            String keyFile = dates.sampleId() + dates.tumorBarcode() + dates.reportDate() + dates.sourceReport();
            if (keySample.equals(keyFile)) {
                LOGGER.warn("Sample is already reported");
                present = true;
            }
        }

        if (!present && addToFile) {
            LOGGER.info("Writing report date to tsv file");
            String stringForFile = sampleId + "\t" + tumorBarcode + "\t" + reportDate + "\t" + reason + "\n";
            writeToTSV(stringForFile, reportDatesTSV);
        }
    }

    private static String currentDate() {
        SimpleDateFormat formatter = new SimpleDateFormat("dd/MM/yyyy");
        return formatter.format(new Date());
    }

    private static void writeToTSV(@NotNull String stringForFile, @NotNull String reportDatesTSV) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(reportDatesTSV, true));
        writer.write(stringForFile);
        writer.close();
    }
}
