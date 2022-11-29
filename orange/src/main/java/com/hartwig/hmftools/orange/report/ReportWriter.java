package com.hartwig.hmftools.orange.report;

import java.io.BufferedWriter;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import com.google.gson.GsonBuilder;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.algo.util.OrangeReportModifier;
import com.hartwig.hmftools.orange.report.chapters.CohortComparisonChapter;
import com.hartwig.hmftools.orange.report.chapters.FrontPageChapter;
import com.hartwig.hmftools.orange.report.chapters.GermlineFindingsChapter;
import com.hartwig.hmftools.orange.report.chapters.ImmunologyChapter;
import com.hartwig.hmftools.orange.report.chapters.QualityControlChapter;
import com.hartwig.hmftools.orange.report.chapters.RNAFindingsChapter;
import com.hartwig.hmftools.orange.report.chapters.ReportChapter;
import com.hartwig.hmftools.orange.report.chapters.SomaticFindingsChapter;
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
import org.jetbrains.annotations.Nullable;

public class ReportWriter {

    private static final Logger LOGGER = LogManager.getLogger(ReportWriter.class);

    private final boolean writeToDisk;
    @Nullable
    private final String outputDir;
    @NotNull
    private final ReportConfig reportConfig;

    ReportWriter(final boolean writeToDisk, @Nullable final String outputDir, @NotNull final ReportConfig reportConfig) {
        this.writeToDisk = writeToDisk;
        this.outputDir = outputDir;
        this.reportConfig = reportConfig;
    }

    public void write(@NotNull OrangeReport report) throws IOException {
        writePdf(report);
        writeJson(report);
    }

    private void writePdf(@NotNull OrangeReport report) throws IOException {
        ReportChapter[] chapters = new ReportChapter[] { new FrontPageChapter(report, reportConfig.reportGermline()),
                new SomaticFindingsChapter(report), new GermlineFindingsChapter(report, reportConfig.reportGermline()),
                new ImmunologyChapter(report), new RNAFindingsChapter(report), new CohortComparisonChapter(report),
                new QualityControlChapter(report) };

        String platinumVersion = report.platinumVersion() != null ? report.platinumVersion() : ReportResources.NOT_AVAILABLE;
        writePdfChapters(report.sampleId(), platinumVersion, chapters);
    }

    private void writeJson(@NotNull OrangeReport report) throws IOException {
        if (writeToDisk && outputDir != null) {
            String outputFilePath = outputDir + File.separator + report.sampleId() + ".orange.json";
            LOGGER.info("Writing JSON report to {} ", outputFilePath);

            OrangeReport reportToWrite = reportConfig.limitJsonOutput() ? OrangeReportModifier.limitAllListsToMaxOne(report) : report;
            String json = new GsonBuilder().serializeNulls().serializeSpecialFloatingPointValues().create().toJson(reportToWrite);
            BufferedWriter writer = new BufferedWriter(new FileWriter(outputFilePath));

            writer.write(json);
            writer.close();
        } else {
            LOGGER.info("Generating in-memory JSON report");
        }
    }

    private void writePdfChapters(@NotNull String sampleId, @NotNull String platinumVersion, @NotNull ReportChapter[] chapters)
            throws IOException {
        Document doc = initializeReport(sampleId);
        PdfDocument pdfDocument = doc.getPdfDocument();

        PageEventHandler pageEventHandler = PageEventHandler.create(sampleId, platinumVersion);
        pdfDocument.addEventHandler(PdfDocumentEvent.START_PAGE, pageEventHandler);

        for (int i = 0; i < chapters.length; i++) {
            ReportChapter chapter = chapters[i];
            pdfDocument.setDefaultPageSize(chapter.pageSize());
            pageEventHandler.chapterTitle(chapter.name());
            pageEventHandler.resetChapterPageCounter();

            if (i > 0) {
                doc.add(new AreaBreak(AreaBreakType.NEXT_PAGE));
            }
            chapter.render(doc);
        }

        pageEventHandler.writeTotalPageCount(doc.getPdfDocument());

        doc.close();
        pdfDocument.close();
    }

    @NotNull
    private Document initializeReport(@NotNull String sampleId) throws IOException {
        PdfWriter writer;
        if (writeToDisk) {
            String outputFilePath = outputDir + File.separator + sampleId + ".orange.pdf";
            LOGGER.info("Writing PDF report to {}", outputFilePath);
            writer = new PdfWriter(outputFilePath);
        } else {
            LOGGER.info("Generating in-memory PDF report");
            writer = new PdfWriter(new ByteArrayOutputStream());
        }

        PdfDocument pdf = new PdfDocument(writer);
        pdf.setDefaultPageSize(PageSize.A4);
        pdf.getDocumentInfo().setTitle(ReportResources.METADATA_TITLE);
        pdf.getDocumentInfo().setAuthor(ReportResources.METADATA_AUTHOR);

        Document document = new Document(pdf);
        document.setMargins(ReportResources.PAGE_MARGIN_TOP,
                ReportResources.PAGE_MARGIN_RIGHT,
                ReportResources.PAGE_MARGIN_BOTTOM,
                ReportResources.PAGE_MARGIN_LEFT);

        return document;
    }
}
