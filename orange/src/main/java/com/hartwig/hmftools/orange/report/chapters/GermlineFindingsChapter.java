package com.hartwig.hmftools.orange.report.chapters;

import java.io.IOException;
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
import com.hartwig.hmftools.orange.algo.QcStatusInterpretation;
import com.hartwig.hmftools.orange.report.DocumentContext;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.tables.DisruptionTable;
import com.hartwig.hmftools.orange.report.tables.GainDeletionTable;
import com.hartwig.hmftools.orange.report.tables.GermlineVariantTable;
import com.hartwig.hmftools.orange.report.tables.PharmacogeneticsTable;
import com.hartwig.hmftools.orange.report.util.Cells;

import org.apache.pdfbox.pdmodel.PDPage;
import org.apache.pdfbox.pdmodel.common.PDRectangle;
import org.jetbrains.annotations.NotNull;

import be.quodlibet.boxable.BaseTable;
import be.quodlibet.boxable.Row;

public class GermlineFindingsChapter implements ReportChapter
{
    private final OrangeRecord mReport;
    private final ReportResources mReportResources;

    public GermlineFindingsChapter(final OrangeRecord report, final ReportResources reportResources)
    {
        mReport = report;
        mReportResources = reportResources;
    }

    @NotNull
    @Override
    public String name()
    {
        return "Germline Findings";
    }

    @NotNull
    @Override
    public PDRectangle pageSize()
    {
        return PDRectangle.A4;
    }

    @Override
    public void render(@NotNull final DocumentContext document) throws IOException
    {
        document.addParagraph(name(), mReportResources.chapterTitleStyle());

        if(QcStatusInterpretation.hasPurpleFail(mReport.purple().fit().qc()))
        {
            document.addQcFailNotice(mReportResources);
            return;
        }

        if(mReport.referenceId() != null)
        {
            addGermlineVariants(document);
            document.addSpacing(10);
            addGermlineAmpDels(document);
            document.addSpacing(10);
            addGermlineBreakends(document);
            document.addSpacing(10);
            addGermlineCNAberrations(document);
            document.addSpacing(10);
            addPharmacogenetics(document);
        }
        else
        {
            document.addParagraph(ReportResources.NOT_AVAILABLE, mReportResources.tableContentStyle());
        }
    }

    private void addGermlineVariants(final DocumentContext document) throws IOException
    {
        List<PurpleDriver> drivers = mReport.purple().germlineDrivers();

        List<PurpleVariant> reportableVariants = mReport.purple().germlineVariants();
        if(drivers != null && reportableVariants != null)
        {
            String titleDrivers = "Small Variants (" + reportableVariants.size() + ")";
            document.addTable(GermlineVariantTable.build(document, titleDrivers, contentWidth(), reportableVariants, mReportResources));
        }
    }

    private void addGermlineAmpDels(final DocumentContext document) throws IOException
    {
        List<PurpleGainDeletion> reportableGermlineFullDels = mReport.purple().germlineGainsDels();
        if(reportableGermlineFullDels != null)
        {
            String title = "Amplifications and Deletions (" + reportableGermlineFullDels.size() + ")";
            document.addTable(GainDeletionTable.build(
                    document, title, contentWidth(), reportableGermlineFullDels, mReportResources, mReport.hasRna()));
        }
    }

    private void addGermlineBreakends(final DocumentContext document) throws IOException
    {
        List<LinxBreakend> germlineBreakends = mReport.linx().germlineBreakends();

        if(germlineBreakends != null)
        {
            String title = "Disruptions (" + germlineBreakends.size() + ")";
            document.addTable(DisruptionTable.build(document, title, contentWidth(), germlineBreakends, mReportResources));
        }
    }

    private void addGermlineCNAberrations(final DocumentContext document) throws IOException
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

            Cells cells = new Cells(mReportResources);
            BaseTable table = document.createTable(contentWidth(), null);
            Row<PDPage> titleRow = table.createRow(15f);
            cells.applyTitleStyle(titleRow.createCell(100, "Chromosomal aberrations (" + count + ")"));
            Row<PDPage> dataRow = table.createRow(12f);
            cells.addContentCell(dataRow, 100, germlineAberrationJoiner.toString());
            document.addTable(table);
        }
    }

    private void addPharmacogenetics(final DocumentContext document) throws IOException
    {
        Set<PeachGenotype> peach = mReport.peach();
        if(peach != null)
        {
            String titlePharmacogenetics = "Pharmacogenetics (" + peach.size() + ")";
            document.addTable(PharmacogeneticsTable.build(document, titlePharmacogenetics, contentWidth(), peach, mReportResources));
        }
    }
}
