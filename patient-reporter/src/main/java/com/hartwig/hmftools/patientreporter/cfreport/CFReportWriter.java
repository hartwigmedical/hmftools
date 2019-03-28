package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.QCFailReport;
import com.hartwig.hmftools.patientreporter.ReportWriter;
import com.hartwig.hmftools.patientreporter.cfreport.chapters.*;
import com.hartwig.hmftools.patientreporter.cfreport.components.Footer;

import com.itextpdf.kernel.events.PdfDocumentEvent;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.kernel.pdf.PdfDocument;
import com.itextpdf.kernel.pdf.PdfWriter;
import com.itextpdf.layout.Document;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

public class CFReportWriter implements ReportWriter {

    private static final Logger LOGGER = LogManager.getLogger(CFReportWriter.class);

    @Override
    public void writeAnalysedPatientReport(@NotNull final AnalysedPatientReport report, @NotNull final String outputFilePath) throws IOException {

        // Initialize report with metadata
        final Document reportDocument = initializeReport(outputFilePath);
        final PdfDocument pdfDocument = reportDocument.getPdfDocument();

        // Setup page event handling (used for automatic generation of header, side panel and footer)
        final PageEventHandler pageEventHandler = new PageEventHandler(report);
        pdfDocument.addEventHandler(PdfDocumentEvent.START_PAGE, pageEventHandler);

        // Add chapters
        new SummaryChapter().render(pageEventHandler, report, reportDocument);
        new TherapyDetailsChapter().render(pageEventHandler, report, reportDocument);
        new ActionableOrDriversChapter().render(pageEventHandler, report, reportDocument);
        new TumorCharacteristicsChapter().render(pageEventHandler, report, reportDocument);
        new CircosChapter().render(pageEventHandler, report, reportDocument);
        new ExplanationChapter().render(pageEventHandler, report, reportDocument);
        new DetailsAndDisclaimerChapter().render(pageEventHandler, report, reportDocument);

        // Update total page count on pages and close document
        pageEventHandler.writeTotalPageCount(pdfDocument);
        reportDocument.close();
        pdfDocument.close();

        LOGGER.info("Created patient report at " + outputFilePath);

    }

    @Override
    public void writeQCFailReport(@NotNull final QCFailReport report, @NotNull final String outputFilePath) throws IOException {
        // Does not need to be implemented
    }

    /**
     * Initialize report document with relevant metadata
     */
    @NotNull
    private static Document initializeReport(@NotNull final String outputFilePath) throws IOException {

        // Prevent overwriting existing files
        if (Files.exists(new File(outputFilePath).toPath())) {
            throw new IOException("Could not write " + outputFilePath + " as it already exists.");
        }

        // Create PDF with metadata
        final PdfDocument pdf = new PdfDocument(new PdfWriter(outputFilePath));
        pdf.setDefaultPageSize(PageSize.A4);
        pdf.getDocumentInfo().setTitle(ReportResources.METADATA_TITLE);
        pdf.getDocumentInfo().setAuthor(ReportResources.METADATA_AUTHOR);

        // Create document
        final Document document = new Document(pdf);
        document.setMargins(ReportResources.PAGE_MARGIN_TOP,
                ReportResources.PAGE_MARGIN_RIGHT,
                ReportResources.PAGE_MARGIN_BOTTOM,
                ReportResources.PAGE_MARGIN_LEFT);

        return document;
    }

}
