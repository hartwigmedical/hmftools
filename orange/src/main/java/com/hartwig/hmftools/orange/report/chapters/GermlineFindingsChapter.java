package com.hartwig.hmftools.orange.report.chapters;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleGermlineAberration;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.datamodel.BreakendEntry;
import com.hartwig.hmftools.orange.report.datamodel.BreakendEntryFactory;
import com.hartwig.hmftools.orange.report.datamodel.VariantEntry;
import com.hartwig.hmftools.orange.report.datamodel.VariantEntryFactory;
import com.hartwig.hmftools.orange.report.tables.DisruptionTable;
import com.hartwig.hmftools.orange.report.tables.GainDeletionTable;
import com.hartwig.hmftools.orange.report.tables.GermlineVariantTable;
import com.hartwig.hmftools.orange.report.tables.PharmacogeneticsTable;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;

public class GermlineFindingsChapter implements ReportChapter
{
    private final OrangeRecord mReport;
    private final ReportResources mReportResources;

    public GermlineFindingsChapter(final OrangeRecord report, final ReportResources reportResources)
    {
        mReport = report;
        mReportResources = reportResources;
    }

    @Override
    public String name()
    {
        return "Germline Findings";
    }

    @Override
    public PageSize pageSize()
    {
        return PageSize.A4;
    }

    @Override
    public void render(final Document document)
    {
        document.add(new Paragraph(name()).addStyle(mReportResources.chapterTitleStyle()));

        if(mReport.referenceId() != null)
        {
            addGermlineVariants(document);
            addGermlineAmpDels(document);
            addGermlineBreakends(document);
            addGermlineCNAberrations(document);
            addPharmacogenetics(document);
        }
        else
        {
            document.add(new Paragraph(ReportResources.NOT_AVAILABLE).addStyle(mReportResources.tableContentStyle()));
        }
    }

    private void addGermlineVariants(final Document document)
    {
        List<PurpleDriver> drivers = mReport.purple().germlineDrivers();

        List<PurpleVariant> reportableVariants = mReport.purple().germlineVariants();
        if(drivers != null && reportableVariants != null)
        {
            List<VariantEntry> reportableEntries = VariantEntryFactory.create(reportableVariants, drivers);
            String titleDrivers = "Small variants (" + reportableEntries.size() + ")";
            document.add(GermlineVariantTable.build(titleDrivers, contentWidth(), reportableEntries, mReportResources));
        }
    }

    private void addGermlineAmpDels(final Document document)
    {
        List<PurpleGainDeletion> reportableGermlineFullDels = mReport.purple().germlineGainsDels();
        if(reportableGermlineFullDels != null)
        {
            String title = "Amplifications and Deletions (" + reportableGermlineFullDels.size() + ")";

            document.add(GainDeletionTable.build(
                    title, contentWidth(), reportableGermlineFullDels, mReportResources, mReport.hasRna()));
        }
    }

    private void addGermlineBreakends(final Document document)
    {
        List<LinxBreakend> reportableGermlineBreakends = mReport.linx().germlineBreakends();

        if(reportableGermlineBreakends != null)
        {
            List<BreakendEntry> reportableBreakends = BreakendEntryFactory.create(reportableGermlineBreakends);

            String title = "Disruptions (" + reportableBreakends.size() + ")";
            document.add(DisruptionTable.build(title, contentWidth(), reportableBreakends, mReportResources));
        }
    }

    private void addGermlineCNAberrations(final Document document)
    {
        Set<PurpleGermlineAberration> germlineAberrations = mReport.purple().fit().qc().germlineAberrations();
        if(!germlineAberrations.isEmpty())
        {
            int count = 0;
            StringJoiner germlineAberrationJoiner = new StringJoiner(", ");
            for(PurpleGermlineAberration germlineAberration : germlineAberrations)
            {
                if(germlineAberration != PurpleGermlineAberration.NONE)
                {
                    count++;
                }
                germlineAberrationJoiner.add(germlineAberration.toString());
            }
            Table table = new Table(UnitValue.createPercentArray(new float[] { 1 })).setWidth(contentWidth());
            table.addCell(new Cells(mReportResources).createContent(germlineAberrationJoiner.toString()));
            document.add(new Tables(mReportResources).createWrapping(table, "Germline CN aberrations (" + count + ")"));
        }
    }

    private void addPharmacogenetics(final Document document)
    {
        Set<PeachGenotype> peach = mReport.peach();
        if(peach != null)
        {
            String titlePharmacogenetics = "Pharmacogenetics (" + peach.size() + ")";
            document.add(PharmacogeneticsTable.build(titlePharmacogenetics, contentWidth(), peach, mReportResources));
        }
    }
}
