package com.hartwig.hmftools.orange.report;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.datamodel.OrangeJson;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.report.chapters.CohortComparisonChapter;
import com.hartwig.hmftools.orange.report.chapters.FrontPageChapter;
import com.hartwig.hmftools.orange.report.chapters.GermlineFindingsChapter;
import com.hartwig.hmftools.orange.report.chapters.ImmunologyChapter;
import com.hartwig.hmftools.orange.report.chapters.QualityControlChapter;
import com.hartwig.hmftools.orange.report.chapters.ReportChapter;
import com.hartwig.hmftools.orange.report.chapters.RnaFindingsChapter;
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

import org.jetbrains.annotations.Nullable;

public class ReportWriter
{
    private final boolean mWriteToDisk;

    @Nullable
    private final String mOutputDir;
    private final PlotPathResolver mPlotPathResolver;
    private final boolean mAddDisclaimer;

    ReportWriter(boolean writeToDisk, @Nullable String outputDir, final PlotPathResolver plotPathResolver, boolean addDisclaimer)
    {
        mWriteToDisk = writeToDisk;
        mOutputDir = outputDir;
        mPlotPathResolver = plotPathResolver;
        mAddDisclaimer = addDisclaimer;
    }

    public void write(final OrangeRecord report) throws IOException
    {
        writePdf(report);
        writeJson(report);
    }

    private void writePdf(final OrangeRecord report) throws IOException
    {
        ReportResources reportResources = ReportResources.create();

        List<ReportChapter> chapters = new ArrayList<>();

        chapters.add(new FrontPageChapter(report, mPlotPathResolver, reportResources));
        chapters.add(new SomaticFindingsChapter(report, mPlotPathResolver, reportResources));

        if(!report.tumorOnlyMode())
        {
            chapters.add(new GermlineFindingsChapter(report, reportResources));
        }

        chapters.add(new ImmunologyChapter(report, reportResources));

        IsofoxRecord isofox = report.isofox();
        if(isofox != null)
        {
            chapters.add(new RnaFindingsChapter(isofox, report.purple(), reportResources));
        }

        if(!report.tumorOnlyMode())
        {
            chapters.add(new CohortComparisonChapter(report, mPlotPathResolver, reportResources));
        }

        chapters.add(new QualityControlChapter(report, mPlotPathResolver, reportResources));

        String pipelineVersion = report.pipelineVersion() != null ? report.pipelineVersion() : ReportResources.NOT_AVAILABLE;
        writePdfChapters(report.sampleId(), pipelineVersion, chapters, reportResources);
    }

    private void writeJson(final OrangeRecord report) throws IOException
    {
        if(mWriteToDisk && mOutputDir != null)
        {
            String basePath = FileWriterUtils.checkAddDirSeparator(mOutputDir);
            String outputFilePath = basePath + report.sampleId() + ".orange.json";
            LOGGER.info("Writing JSON report to {} ", outputFilePath);

            OrangeJson.getInstance().write(report, outputFilePath);
        }
        else
        {
            LOGGER.info("Generating in-memory JSON report");
        }
    }

    private void writePdfChapters(
            final String sampleId, final String pipelineVersion, final List<ReportChapter> chapters,
            final ReportResources reportResources) throws IOException
    {
        Document doc = initializeReport(sampleId);
        PdfDocument pdfDocument = doc.getPdfDocument();

        PageEventHandler pageEventHandler = PageEventHandler.create(sampleId, pipelineVersion, reportResources, mAddDisclaimer);
        pdfDocument.addEventHandler(PdfDocumentEvent.START_PAGE, pageEventHandler);

        for(int i = 0; i < chapters.size(); i++)
        {
            ReportChapter chapter = chapters.get(i);
            pdfDocument.setDefaultPageSize(chapter.pageSize());
            pageEventHandler.chapterTitle(chapter.name());
            pageEventHandler.resetChapterPageCounter();

            if(i > 0)
            {
                doc.add(new AreaBreak(AreaBreakType.NEXT_PAGE));
            }
            chapter.render(doc);
        }

        pageEventHandler.writeFooters(doc.getPdfDocument());

        doc.close();
        pdfDocument.close();
    }

    private Document initializeReport(final String sampleId) throws IOException
    {
        PdfWriter writer;
        WriterProperties properties = new WriterProperties()
                .setFullCompressionMode(true)
                .setCompressionLevel(CompressionConstants.BEST_COMPRESSION)
                .useSmartMode();
        if(mWriteToDisk)
        {
            String basePath = FileWriterUtils.checkAddDirSeparator(mOutputDir);
            String outputFilePath = basePath + sampleId + ".orange.pdf";
            LOGGER.info("Writing PDF report to {}", outputFilePath);
            writer = new PdfWriter(outputFilePath, properties);
        }
        else
        {
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
