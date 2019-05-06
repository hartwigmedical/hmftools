package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.QCFailReport;
import com.hartwig.hmftools.patientreporter.ReportWriter;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.cfreport.chapters.*;
import com.itextpdf.kernel.events.PdfDocumentEvent;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.kernel.pdf.PdfDocument;
import com.itextpdf.kernel.pdf.PdfWriter;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.AreaBreak;
import com.itextpdf.layout.property.AreaBreakType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

public class CFReportWriter implements ReportWriter {

    private static final Logger LOGGER = LogManager.getLogger(CFReportWriter.class);

    private boolean writeToFile = true;

    @Override
    public void writeAnalysedPatientReport(@NotNull final AnalysedPatientReport report, @NotNull final String outputFilePath)
            throws IOException {
        writeReport(report.sampleReport(),
                new ReportChapter[] { new SummaryChapter(report), new TherapyDetailsChapter(report), new ActionableOrDriversChapter(report),
                        new TumorCharacteristicsChapter(report), new CircosChapter(report), new ExplanationChapter(),
                        new DetailsAndDisclaimerChapter(report) },
                outputFilePath);
    }

    @Override
    public void writeQCFailReport(@NotNull final QCFailReport report, @NotNull final String outputFilePath) throws IOException {
        writeReport(report.sampleReport(), new ReportChapter[] { new QCFailChapter(report) }, outputFilePath);
    }

    public void setWriteToFile(boolean writeToFile) {
        this.writeToFile = writeToFile;
    }

    private void writeReport(@NotNull final SampleReport sampleReport, @NotNull ReportChapter[] chapters, @NotNull String outputFilePath)
            throws IOException {
        // Initialize report with metadata
        final Document doc = initializeReport(outputFilePath);
        final PdfDocument pdfDocument = doc.getPdfDocument();

        // Setup page event handling (used for automatic generation of header, side panel and footer)
        PageEventHandler pageEventHandler = new PageEventHandler(sampleReport);
        pdfDocument.addEventHandler(PdfDocumentEvent.START_PAGE, pageEventHandler);

        // Write chapters to report
        for (int i = 0; i < chapters.length; i++) {
            ReportChapter chapter = chapters[i];

            // Reconfigure event handler for upcoming pages
            pageEventHandler.setChapterTitle(chapter.getName());
            pageEventHandler.setPageNumberPrefix(chapter.getPageNumberPrefix());
            pageEventHandler.resetChapterPageCounter();
            pageEventHandler.setSidebarType(!chapter.isFullWidth(), chapter.hasCompleteSidebar());

            // Add page break for all except the first chapters
            if (i > 0) {
                doc.add(new AreaBreak(AreaBreakType.NEXT_PAGE));
            }
            chapter.render(doc);

        }

        // Update chapter titles and page numbers
        pageEventHandler.writeDynamicTextParts(doc.getPdfDocument());

        // Close document
        doc.close();
        pdfDocument.close();

        // Log info
        if (writeToFile) {
            LOGGER.info("Created patient report at " + outputFilePath);
        } else {
            LOGGER.info("Created patient report");
        }
    }

    /**
     * Initialize report document with relevant metadata
     */
    @NotNull
    private Document initializeReport(@NotNull final String outputFilePath) throws IOException {
        PdfWriter writer;
        if (writeToFile) {
            // Prevent overwriting existing files
            if (Files.exists(new File(outputFilePath).toPath())) {
                throw new IOException("Could not write " + outputFilePath + " as it already exists.");
            }

            // Create file writer
            writer = new PdfWriter(outputFilePath);

        } else {
            // Create in-memory writer
            writer = new PdfWriter(new ByteArrayOutputStream());
        }

        // Create PDF with metadata
        final PdfDocument pdf = new PdfDocument(writer);
        pdf.setDefaultPageSize(PageSize.A4);
        pdf.getDocumentInfo().setTitle(ReportResources.METADATA_TITLE);
        pdf.getDocumentInfo().setAuthor(ReportResources.METADATA_AUTHOR);

        // Create document
        Document document = new Document(pdf);
        document.setMargins(ReportResources.PAGE_MARGIN_TOP,
                ReportResources.PAGE_MARGIN_RIGHT,
                ReportResources.PAGE_MARGIN_BOTTOM,
                ReportResources.PAGE_MARGIN_LEFT);

        return document;
    }
}
