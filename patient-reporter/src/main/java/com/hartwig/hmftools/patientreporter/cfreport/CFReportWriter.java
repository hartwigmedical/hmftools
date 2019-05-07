package com.hartwig.hmftools.patientreporter.cfreport;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.QCFailReport;
import com.hartwig.hmftools.patientreporter.ReportWriter;
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

    private final boolean writeToFile;

    @NotNull
    public static CFReportWriter createProductionReportWriter() {
        return new CFReportWriter(true);
    }

    @VisibleForTesting
    CFReportWriter(final boolean writeToFile) {
        this.writeToFile = writeToFile;
    }

    @Override
    public void writeAnalysedPatientReport(@NotNull final AnalysedPatientReport report, @NotNull final String outputFilePath)
            throws IOException {
        writeReport(report,
                new ReportChapter[] { new SummaryChapter(report), new TherapyDetailsChapter(report), new ActionableOrDriversChapter(report),
                        new TumorCharacteristicsChapter(report), new CircosChapter(report), new ExplanationChapter(),
                        new DetailsAndDisclaimerChapter(report) },
                outputFilePath);
    }

    @Override
    public void writeQCFailReport(@NotNull final QCFailReport report, @NotNull final String outputFilePath) throws IOException {
        writeReport(report, new ReportChapter[] { new QCFailChapter(report) }, outputFilePath);
    }

    private void writeReport(@NotNull final PatientReport patientReport, @NotNull ReportChapter[] chapters, @NotNull String outputFilePath)
            throws IOException {
        final Document doc = initializeReport(outputFilePath, writeToFile);
        final PdfDocument pdfDocument = doc.getPdfDocument();

        final PageEventHandler pageEventHandler = new PageEventHandler(patientReport);
        pdfDocument.addEventHandler(PdfDocumentEvent.START_PAGE, pageEventHandler);

        for (int i = 0; i < chapters.length; i++) {
            ReportChapter chapter = chapters[i];

            pageEventHandler.setChapterTitle(chapter.name());
            pageEventHandler.setPageNumberPrefix(chapter.pageNumberPrefix());
            pageEventHandler.resetChapterPageCounter();
            pageEventHandler.setSidebarType(!chapter.isFullWidth(), chapter.hasCompleteSidebar());

            if (i > 0) {
                doc.add(new AreaBreak(AreaBreakType.NEXT_PAGE));
            }
            chapter.render(doc);
        }

        pageEventHandler.writeDynamicTextParts(doc.getPdfDocument());

        doc.close();
        pdfDocument.close();

        if (writeToFile) {
            LOGGER.info("Created patient report at " + outputFilePath);
        } else {
            LOGGER.info("Created patient report");
        }
    }

    @NotNull
    private static Document initializeReport(@NotNull String outputFilePath, boolean writeToFile) throws IOException {
        PdfWriter writer;
        if (writeToFile) {
            if (Files.exists(new File(outputFilePath).toPath())) {
                throw new IOException("Could not write " + outputFilePath + " as it already exists.");
            }

            writer = new PdfWriter(outputFilePath);
        } else {
            // Write output to output stream where it is effectively ignored.
            writer = new PdfWriter(new ByteArrayOutputStream());
        }

        final PdfDocument pdf = new PdfDocument(writer);
        pdf.setDefaultPageSize(PageSize.A4);
        pdf.getDocumentInfo().setTitle(ReportResources.METADATA_TITLE);
        pdf.getDocumentInfo().setAuthor(ReportResources.METADATA_AUTHOR);

        final Document document = new Document(pdf);
        document.setMargins(ReportResources.PAGE_MARGIN_TOP,
                ReportResources.PAGE_MARGIN_RIGHT,
                ReportResources.PAGE_MARGIN_BOTTOM,
                ReportResources.PAGE_MARGIN_LEFT);

        return document;
    }
}
