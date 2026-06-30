package com.hartwig.hmftools.orange.report;

import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;

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
import com.hartwig.hmftools.orange.report.pdfdata.CuppaChapterData;
import com.hartwig.hmftools.orange.report.pdfdata.CuppaChapterDataFactory;
import com.hartwig.hmftools.orange.report.pdfdata.FrontPageData;
import com.hartwig.hmftools.orange.report.pdfdata.FrontPageDataFactory;
import com.hartwig.hmftools.orange.report.pdfdata.GermlineFindingsData;
import com.hartwig.hmftools.orange.report.pdfdata.GermlineFindingsDataFactory;
import com.hartwig.hmftools.orange.report.pdfdata.ImmunologyData;
import com.hartwig.hmftools.orange.report.pdfdata.ImmunologyDataFactory;
import com.hartwig.hmftools.orange.report.pdfdata.PurplePlotsData;
import com.hartwig.hmftools.orange.report.pdfdata.PurplePlotsDataFactory;
import com.hartwig.hmftools.orange.report.pdfdata.QualityControlData;
import com.hartwig.hmftools.orange.report.pdfdata.QualityControlDataFactory;
import com.hartwig.hmftools.orange.report.pdfdata.RnaFindingsData;
import com.hartwig.hmftools.orange.report.pdfdata.RnaFindingsDataFactory;
import com.hartwig.hmftools.orange.report.pdfdata.SomaticFindingsData;
import com.hartwig.hmftools.orange.report.pdfdata.SomaticFindingsDataFactory;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.DynamicReports;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;
import net.sf.jasperreports.engine.JREmptyDataSource;

public class OrangeReport
{
    private final OrangeConfig mConfig;
    private final PlotPathResolver mPlotPathResolver;

    public OrangeReport(final OrangeConfig config, final PlotPathResolver plotPathResolver)
    {
        mConfig = config;
        mPlotPathResolver = plotPathResolver;
    }

    public void writeReport(final OrangeRecord orangeData, OutputStream output)
    {
        String displaySampleId = mConfig != null ? mConfig.DisplaySampleId : orangeData.sampleId();
        List<ReportChapter> chapters = buildChapters(orangeData);

        List<JasperReportBuilder> chapterReports = new ArrayList<>();
        for(ReportChapter chapter : chapters)
        {
            JasperReportBuilder chapterReport = chapter.buildReport();
            //            OrangeTemplateFactory.applyHeaderLayout(chapterReport, displaySampleId);
            chapterReports.add(chapterReport);
        }
        VerticalListBuilder combinedContent = DynamicReports.cmp.verticalList();
        for(int i = 0; i < chapters.size(); i++)
        {
            if(i > 0)
            {
                combinedContent.add(DynamicReports.cmp.pageBreak());
            }
            combinedContent.add(DynamicReports.cmp.subreport(chapters.get(i).buildReport()));
        }
        try
        {
            DynamicReports.report()
                    // Define global Page Header & Footer here so they repeat on ALL chapter pages
                    .pageHeader(OrangeTemplateFactory.header(displaySampleId))
                    .pageFooter(OrangeTemplateFactory.footer())
                    .setSummaryWithPageHeaderAndFooter(true)
                    // Add the combined chapters into the summary band
                    .summary(combinedContent)

                    // Master report needs at least a dummy 1-row data source to trigger rendering
                    .setDataSource(new JREmptyDataSource())
                    .toPdf(output);
        }
        catch(Exception e)
        {
            throw new RuntimeException(e);
        }
    }

    private List<ReportChapter> buildChapters(final OrangeRecord report)
    {
        List<ReportChapter> chapters = new ArrayList<>();

        FrontPageData frontPageData = FrontPageDataFactory.build(report, mConfig, mPlotPathResolver);
        chapters.add(new FrontPageChapter(frontPageData, null));

        SomaticFindingsData somaticData = SomaticFindingsDataFactory.build(report, mConfig, mPlotPathResolver);
        chapters.add(new SomaticFindingsChapter(somaticData, null));

        if(!report.tumorOnlyMode())
        {
            GermlineFindingsData germlineData = GermlineFindingsDataFactory.build(report);
            chapters.add(new GermlineFindingsChapter(germlineData, null));
        }

        ImmunologyData immunologyData = ImmunologyDataFactory.build(report);
        chapters.add(new ImmunologyChapter(immunologyData, null));

        IsofoxRecord isofox = report.isofox();
        if(isofox != null)
        {
            RnaFindingsData rnaData = RnaFindingsDataFactory.build(isofox);
            chapters.add(new RnaFindingsChapter(rnaData, null));
        }

        if(!report.tumorOnlyMode())
        {
            CuppaChapterData cuppaData = CuppaChapterDataFactory.build(report, mPlotPathResolver);
            chapters.add(new CuppaChapter(cuppaData, null));
        }

        if(report.plots().qSeePlot() != null)
        {
            QualityControlData qcData = QualityControlDataFactory.build(report, mPlotPathResolver);
            chapters.add(new QualityControlChapter(qcData, null));
        }

        PurplePlotsData purplePlotsData = PurplePlotsDataFactory.build(report, mPlotPathResolver);
        chapters.add(new PurplePlotsChapter(purplePlotsData, null));

        return chapters;
    }

}
