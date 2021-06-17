package com.hartwig.hmftools.orange.report;

import com.hartwig.hmftools.orange.algo.OrangeReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ReportWriter {

    private static final Logger LOGGER = LogManager.getLogger(ReportWriter.class);

    @NotNull
    private final String outputDir;

    public ReportWriter(@NotNull final String outputDir) {
        this.outputDir = outputDir;
    }

    public void write(@NotNull OrangeReport report) {
        LOGGER.info("Writing {}", report);
//        Document doc = initializeReport(outputFilePath, writeToFile);
//        PdfDocument pdfDocument = doc.getPdfDocument();
//
//        PageEventHandler pageEventHandler = new PageEventHandler(patientReport);
//        pdfDocument.addEventHandler(PdfDocumentEvent.START_PAGE, pageEventHandler);
//
//        for (int i = 0; i < chapters.length; i++) {
//            ReportChapter chapter = chapters[i];
//
//            pageEventHandler.chapterTitle(chapter.name());
//            pageEventHandler.resetChapterPageCounter();
//            pageEventHandler.sidebarType(!chapter.isFullWidth(), chapter.hasCompleteSidebar());
//
//            if (i > 0) {
//                doc.add(new AreaBreak(AreaBreakType.NEXT_PAGE));
//            }
//            chapter.render(doc);
//        }
//
//        pageEventHandler.writeDynamicTextParts(doc.getPdfDocument());
//
//        doc.close();
//        pdfDocument.close();
//
//        if (writeToFile) {
//            LOGGER.info("Created patient report at {}", outputFilePath);
//        } else {
//            LOGGER.info("Successfully generated in-memory patient report");
//        }
//    }
//
//    @NotNull
//    private static Document initializeReport(@NotNull String outputFilePath, boolean writeToFile) throws IOException {
//        PdfWriter writer;
//        if (writeToFile) {
//            if (Files.exists(new File(outputFilePath).toPath())) {
//                throw new IOException("Could not write " + outputFilePath + " as it already exists.");
//            }
//
//            writer = new PdfWriter(outputFilePath);
//        } else {
//            // Write output to output stream where it is effectively ignored.
//            writer = new PdfWriter(new ByteArrayOutputStream());
//        }
//
//        PdfDocument pdf = new PdfDocument(writer);
//        pdf.setDefaultPageSize(PageSize.A4);
//        pdf.getDocumentInfo().setTitle(ReportResources.METADATA_TITLE);
//        pdf.getDocumentInfo().setAuthor(ReportResources.METADATA_AUTHOR);
//
//        Document document = new Document(pdf);
//        document.setMargins(ReportResources.PAGE_MARGIN_TOP,
//                ReportResources.PAGE_MARGIN_RIGHT,
//                ReportResources.PAGE_MARGIN_BOTTOM,
//                ReportResources.PAGE_MARGIN_LEFT);
//
//        return document;
    }
}
