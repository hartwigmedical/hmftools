package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.QCFailReport;
import com.hartwig.hmftools.patientreporter.ReportWriter;

import com.lowagie.text.*;
import com.lowagie.text.pdf.PdfWriter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import java.io.*;
import java.nio.file.Files;

public class CFReportWriter implements ReportWriter {

    private static final Logger LOGGER = LogManager.getLogger(CFReportWriter.class);

    @Override
    public void writeAnalysedPatientReport(@NotNull final AnalysedPatientReport report, @NotNull final String outputFilePath) {
        // TODO!
    }

    @Override
    public void writeQCFailReport(@NotNull final QCFailReport report, @NotNull final String outputFilePath) {
        // Does not need to be implemented
    }

    /**
     * Temporary main method for quick running while report is in development
     * @TODO: Remove this
     * @deprecated
     */
    public static void main(String[] args) {

        String outputFilePath = "/Users/Wilco/hmf/tmp/temp_" + String.valueOf(System.currentTimeMillis()) + ".pdf";

        // Initialize report and writer
        final Document report = initializeReport();
        final PdfWriter pdfWriter;
        try {
            pdfWriter = initializePdfWriter(report, outputFilePath);
        } catch (Exception e) {
            e.printStackTrace();
            return;
        }

        // Add content to report
        report.open();
        try {
            report.add(new Phrase("Hello world"));
        } catch (DocumentException e) {
            e.printStackTrace();
            return;
        }

        // Close report (gets written to file by closing)
        report.close();
        pdfWriter.close();

    }

    /**
     * Initialize report with relevant metadata
     * @return
     */
    @NotNull
    private static Document initializeReport() {

        //
        Document document = new Document(PageSize.A4,
                ReportConfiguration.PAGE_MARGIN_LEFT,
                ReportConfiguration.PAGE_MARGIN_RIGHT,
                ReportConfiguration.PAGE_MARGIN_TOP,
                ReportConfiguration.PAGE_MARGIN_BOTTOM);

        // Set document metadata
        document.addTitle(ReportConfiguration.METADATA_TITLE);
        document.addAuthor(ReportConfiguration.METADATA_AUTHOR);

        return document;
    }

    /**
     * Initialize a PdfWriter for the given Document
     *
     * @param report            report document
     * @param outputFilePath    file where the document will be written when the writer is closed
     * @return
     * @throws IOException
     */
    private static PdfWriter initializePdfWriter(@NotNull Document report, @NotNull final String outputFilePath) throws IOException {

        if (Files.exists(new File(outputFilePath).toPath())) {
            throw new IOException("Could not write " + outputFilePath + " as it already exists.");
        }

        try {

            final PdfWriter pdfWriter = PdfWriter.getInstance(report, new FileOutputStream(outputFilePath));
            pdfWriter.setPageEvent(new PageEvent());

            return pdfWriter;

        } catch (DocumentException e) {
            throw new IOException("Could not generate pdf due to iText Document exception: " + e.getMessage());
        } catch (FileNotFoundException e) {
            throw new IOException("Could not generate pdf because the file " + outputFilePath + " could not be created");
        }

    }

 }
