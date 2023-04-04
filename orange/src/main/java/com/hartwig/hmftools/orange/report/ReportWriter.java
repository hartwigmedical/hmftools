package com.hartwig.hmftools.orange.report;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.datamodel.OrangeJson;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
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
import com.itextpdf.kernel.pdf.CompressionConstants;
import com.itextpdf.kernel.pdf.PdfDocument;
import com.itextpdf.kernel.pdf.PdfWriter;
import com.itextpdf.kernel.pdf.WriterProperties;
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
    private final PlotPathResolver plotPathResolver;
    private final boolean addDisclaimer;

    ReportWriter(boolean writeToDisk, @Nullable String outputDir, @NotNull PlotPathResolver plotPathResolver, boolean addDisclaimer) {
        this.writeToDisk = writeToDisk;
        this.outputDir = outputDir;
        this.plotPathResolver = plotPathResolver;
        this.addDisclaimer = addDisclaimer;
    }

    public void write(@NotNull OrangeRecord report) throws IOException {
        writePdf(report);
        writeJson(report);
    }

    private void writePdf(@NotNull OrangeRecord report) throws IOException {
        ReportResources reportResources = ReportResources.create();
        ReportChapter[] chapters = new ReportChapter[] { new FrontPageChapter(report, plotPathResolver, reportResources),
                new SomaticFindingsChapter(report, plotPathResolver, reportResources), new GermlineFindingsChapter(report, reportResources),
                new ImmunologyChapter(report, reportResources), new RNAFindingsChapter(report, reportResources),
                new CohortComparisonChapter(report, plotPathResolver, reportResources),
                new QualityControlChapter(report, plotPathResolver, reportResources) };

        String platinumVersion = report.platinumVersion() != null ? report.platinumVersion() : ReportResources.NOT_AVAILABLE;
        writePdfChapters(report.sampleId(), platinumVersion, chapters, reportResources);
    }

    private void writeJson(@NotNull OrangeRecord report) throws IOException {
        if (writeToDisk && outputDir != null) {
            String outputFilePath = outputDir + File.separator + report.sampleId() + ".orange.json";
            LOGGER.info("Writing JSON report to {} ", outputFilePath);

            OrangeJson.getInstance().write(report, outputFilePath);
        } else {
            LOGGER.info("Generating in-memory JSON report");
        }
    }

    private void writePdfChapters(@NotNull String sampleId, @NotNull String platinumVersion, @NotNull ReportChapter[] chapters,
            @NotNull ReportResources reportResources) throws IOException {
        Document doc = initializeReport(sampleId);
        PdfDocument pdfDocument = doc.getPdfDocument();

        PageEventHandler pageEventHandler = PageEventHandler.create(sampleId, platinumVersion, reportResources, addDisclaimer);
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
        WriterProperties properties = new WriterProperties()
                .setFullCompressionMode(true)
                .setCompressionLevel(CompressionConstants.BEST_COMPRESSION)
                .useSmartMode();
        if (writeToDisk) {
            String outputFilePath = outputDir + File.separator + sampleId + ".orange.pdf";
            LOGGER.info("Writing PDF report to {}", outputFilePath);
            writer = new PdfWriter(outputFilePath, properties);
        } else {
            LOGGER.info("Generating in-memory PDF report");
            writer = new PdfWriter(new ByteArrayOutputStream(), properties);
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
