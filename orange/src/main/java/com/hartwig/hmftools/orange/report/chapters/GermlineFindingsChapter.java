package com.hartwig.hmftools.orange.report.chapters;

import java.util.StringJoiner;

import com.hartwig.hmftools.datamodel.purple.PurpleGermlineAberration;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.pdfdata.GermlineFindingsData;
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

import org.jetbrains.annotations.NotNull;

public class GermlineFindingsChapter implements ReportChapter
{
    private final GermlineFindingsData mData;
    private final ReportResources mReportResources;

    public GermlineFindingsChapter(final GermlineFindingsData data, final ReportResources reportResources)
    {
        mData = data;
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
    public PageSize pageSize()
    {
        return PageSize.A4;
    }

    @Override
    public void render(@NotNull final Document document)
    {
        document.add(new Paragraph(name()).addStyle(mReportResources.chapterTitleStyle()));

        if(mData.hasPurpleFail)
        {
            mReportResources.addQcFailNotice(document);
            return;
        }

        if(mData.hasReferenceId)
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
        if(mData.germlineDrivers != null && mData.germlineVariants != null)
        {
            String titleDrivers = "Small Variants (" + mData.germlineVariants.size() + ")";
            document.add(GermlineVariantTable.build(titleDrivers, contentWidth(), mData.germlineVariants, mReportResources));
        }
    }

    private void addGermlineAmpDels(final Document document)
    {
        if(mData.germlineGainsDels != null)
        {
            String title = "Amplifications and Deletions (" + mData.germlineGainsDels.size() + ")";

            document.add(GainDeletionTable.build(
                    title, contentWidth(), mData.germlineGainsDels, mReportResources, mData.hasRna));
        }
    }

    private void addGermlineBreakends(final Document document)
    {
        if(mData.germlineBreakends != null)
        {
            String title = "Disruptions (" + mData.germlineBreakends.size() + ")";
            document.add(DisruptionTable.build(title, contentWidth(), mData.germlineBreakends, mReportResources));
        }
    }

    private void addGermlineCNAberrations(final Document document)
    {
        if(!mData.germlineAberrations.isEmpty())
        {
            int count = 0;
            StringJoiner germlineAberrationJoiner = new StringJoiner(", ");
            for(PurpleGermlineAberration germlineAberration : mData.germlineAberrations)
            {
                if(germlineAberration != PurpleGermlineAberration.NONE)
                {
                    count++;
                }
                germlineAberrationJoiner.add(germlineAberration.toString());
            }
            Table table = new Table(UnitValue.createPercentArray(new float[] { 1 })).setWidth(contentWidth());
            table.addCell(new Cells(mReportResources).createContent(germlineAberrationJoiner.toString()));
            document.add(new Tables(mReportResources).createWrapping(table, "Chromosomal aberrations (" + count + ")"));
        }
    }

    private void addPharmacogenetics(final Document document)
    {
        if(mData.peach != null)
        {
            String titlePharmacogenetics = "Pharmacogenetics (" + mData.peach.size() + ")";
            document.add(PharmacogeneticsTable.build(titlePharmacogenetics, contentWidth(), mData.peach, mReportResources));
        }
    }
}
