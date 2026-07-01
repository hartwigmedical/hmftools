package com.hartwig.hmftools.orange.report;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.datamodel.OrangeJson;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.report.chapters.CuppaChapter;
import com.hartwig.hmftools.orange.report.chapters.FrontPageChapter;
import com.hartwig.hmftools.orange.report.chapters.GermlineFindingsChapter;
import com.hartwig.hmftools.orange.report.chapters.ImmunologyChapter;
import com.hartwig.hmftools.orange.report.chapters.PurplePlotsChapter;
import com.hartwig.hmftools.orange.report.chapters.QualityControlChapter;
import com.hartwig.hmftools.orange.report.chapters.ReportChapter;
import com.hartwig.hmftools.orange.report.chapters.RnaFindingsChapter;
import com.hartwig.hmftools.orange.report.chapters.SomaticFindingsChapter;

import org.apache.pdfbox.pdmodel.PDDocument;
import org.apache.pdfbox.pdmodel.PDDocumentInformation;
import org.apache.pdfbox.pdmodel.common.PDRectangle;
import org.jetbrains.annotations.Nullable;

public class ReportWriter
{
    private final boolean mWriteToDisk;

    @Nullable
    private final OrangeConfig mConfig;
    private final String mOutputDir;
    private final String mOutputId;

    private final PlotPathResolver mPlotPathResolver;

    public ReportWriter(boolean writeToDisk, @Nullable final OrangeConfig config, final PlotPathResolver plotPathResolver)
    {
        mWriteToDisk = writeToDisk;
        mConfig = config;
        mOutputDir = config != null ? config.OutputDir : null;
        mOutputId = config != null ? config.OutputId : null;
        mPlotPathResolver = plotPathResolver;
    }

    public void write(final OrangeRecord report) throws IOException
    {
        writePdf(report);
        writeJson(report);
    }

    private void writePdf(final OrangeRecord report) throws IOException
    {
        try(PDDocument pdfDocument = new PDDocument())
        {
            ReportResources reportResources = ReportResources.create(pdfDocument);

            List<ReportChapter> chapters = new ArrayList<>();

            chapters.add(new FrontPageChapter(mConfig, report, mPlotPathResolver, reportResources));
            chapters.add(new SomaticFindingsChapter(mConfig, report, mPlotPathResolver, reportResources));

            if(!report.tumorOnlyMode())
            {
                chapters.add(new GermlineFindingsChapter(report, reportResources));
            }

            chapters.add(new ImmunologyChapter(report, reportResources));

            IsofoxRecord isofox = report.isofox();
            if(isofox != null)
            {
                chapters.add(new RnaFindingsChapter(isofox, reportResources));
            }

            if(!report.tumorOnlyMode())
            {
                chapters.add(new CuppaChapter(report, mPlotPathResolver, reportResources));
            }

            if(report.plots().qSeePlot() != null)
            {
                chapters.add(new QualityControlChapter(report, mPlotPathResolver, reportResources));
            }

            chapters.add(new PurplePlotsChapter(report, mPlotPathResolver, reportResources));

            writePdfChapters(pdfDocument, report.sampleId(), chapters, reportResources);
        }
    }

    private String formOutputFile(final String sampleId, final String fileId)
    {
        String filename = mOutputDir + sampleId + ".orange.";

        if(mOutputId != null)
        {
            filename += mOutputId + ".";
        }

        return filename + fileId;
    }

    private void writeJson(final OrangeRecord report) throws IOException
    {
        if(mWriteToDisk && mOutputDir != null)
        {
            String outputFilename = formOutputFile(report.sampleId(), "json");
            LOGGER.info("writing JSON report to {} ", outputFilename);

            OrangeJson.getInstance().write(report, outputFilename);
        }
        else
        {
            LOGGER.info("generating in-memory JSON report");
        }
    }

    private void writePdfChapters(
            final PDDocument pdfDocument, final String sampleId, final List<ReportChapter> chapters,
            final ReportResources reportResources) throws IOException
    {
        // Set metadata
        PDDocumentInformation info = pdfDocument.getDocumentInformation();
        info.setTitle(ReportResources.METADATA_TITLE);
        info.setAuthor(ReportResources.METADATA_AUTHOR);

        // Create document context
        DocumentContext docCtx = new DocumentContext(pdfDocument, PDRectangle.A4,
                ReportResources.PAGE_MARGIN_TOP, ReportResources.PAGE_MARGIN_BOTTOM,
                ReportResources.PAGE_MARGIN_LEFT, ReportResources.PAGE_MARGIN_RIGHT);

        PageEventHandler pageEventHandler = PageEventHandler.create(
                mConfig != null ? mConfig.DisplaySampleId : sampleId, reportResources,
                mConfig != null ? mConfig.AddDisclaimer : false, pdfDocument);

        docCtx.setPageEventHandler(pageEventHandler);

        for(int i = 0; i < chapters.size(); i++)
        {
            ReportChapter chapter = chapters.get(i);
            docCtx.setPageSize(chapter.pageSize());
            pageEventHandler.chapterTitle(chapter.name());
            pageEventHandler.resetChapterPageCounter();

            // Start a new page for each chapter
            docCtx.newPage(chapter.pageSize());
            chapter.render(docCtx);
        }

        // Second pass: write footers with total page count
        docCtx.writeFooters();

        // Save
        if(mWriteToDisk && mOutputDir != null)
        {
            String outputFilename = formOutputFile(sampleId, "pdf");
            LOGGER.info("writing PDF report to {}", outputFilename);
            pdfDocument.save(outputFilename);
        }
        else
        {
            LOGGER.info("generating in-memory PDF report");
            // For in-memory, we still render the document but don't save to file
        }
    }
}
