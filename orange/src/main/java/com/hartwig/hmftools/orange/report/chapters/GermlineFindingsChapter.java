package com.hartwig.hmftools.orange.report.chapters;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.hartwig.hmftools.datamodel.purple.PurpleGermlineAberration;
import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.pdfdata.GermlineFindingsData;
import com.hartwig.hmftools.orange.report.tables.DisruptionTable;
import com.hartwig.hmftools.orange.report.tables.GainDeletionTable;
import com.hartwig.hmftools.orange.report.tables.GermlineVariantTable;
import com.hartwig.hmftools.orange.report.tables.PharmacogeneticsTable;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;
import net.sf.dynamicreports.report.constant.PageOrientation;
import net.sf.dynamicreports.report.constant.PageType;
import net.sf.jasperreports.engine.JREmptyDataSource;

public class GermlineFindingsChapter implements ReportChapter
{
    private final GermlineFindingsData mData;

    public GermlineFindingsChapter(final GermlineFindingsData data, final Object unused)
    {
        mData = data;
    }

    @Override
    public String name()
    {
        return "Germline Findings";
    }

    @Override
    public boolean isLandscape()
    {
        return false;
    }

    @Override
    public JasperReportBuilder buildReport()
    {
        JasperReportBuilder report = report().setPageFormat(PageType.A4, PageOrientation.PORTRAIT);
        VerticalListBuilder content = cmp.verticalList();
        content.add(cmp.text(name()).setStyle(OrangeFonts.CHAPTER_TITLE_STYLE));

        if(mData.hasPurpleFail)
        {
            content.add(cmp.text(ReportResources.NOT_AVAILABLE).setStyle(OrangeFonts.TABLE_CONTENT_STYLE));
            return report.summary(content);
        }

        if(!mData.hasReferenceId)
        {
            content.add(cmp.text(ReportResources.NOT_AVAILABLE).setStyle(OrangeFonts.TABLE_CONTENT_STYLE));
            return report.summary(content);
        }

        content.add(cmp.subreport(buildGermlineVariants()));
        content.add(cmp.subreport(buildGermlineAmpDels()));
        content.add(cmp.subreport(buildGermlineBreakends()));

        if(!mData.germlineAberrations.isEmpty())
        {
            content.add(cmp.subreport(buildChromosomalAberrations()));
        }

        if(mData.peach != null)
        {
            content.add(cmp.subreport(buildPharmacogenetics()));
        }

        return report.summary(content);
    }

    private JasperReportBuilder buildGermlineVariants()
    {
        List<com.hartwig.hmftools.datamodel.purple.PurpleVariant> variants =
                mData.germlineVariants != null ? mData.germlineVariants : List.of();
        String title = "Small Variants (" + variants.size() + ")";
        return GermlineVariantTable.build(title, variants, null);
    }

    private JasperReportBuilder buildGermlineAmpDels()
    {
        List<com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion> gainsDels =
                mData.germlineGainsDels != null ? mData.germlineGainsDels : List.of();
        String title = "Amplifications and Deletions (" + gainsDels.size() + ")";
        return GainDeletionTable.build(title, gainsDels, null, mData.hasRna);
    }

    private JasperReportBuilder buildGermlineBreakends()
    {
        List<com.hartwig.hmftools.datamodel.linx.LinxBreakend> breakends =
                mData.germlineBreakends != null ? mData.germlineBreakends : List.of();
        String title = "Disruptions (" + breakends.size() + ")";
        return DisruptionTable.build(title, breakends, null);
    }

    private JasperReportBuilder buildChromosomalAberrations()
    {
        Set<PurpleGermlineAberration> aberrations = mData.germlineAberrations;
        int count = 0;
        StringJoiner joiner = new StringJoiner(", ");
        for(PurpleGermlineAberration aberration : aberrations)
        {
            if(aberration != PurpleGermlineAberration.NONE)
            {
                count++;
            }
            joiner.add(aberration.toString());
        }
        String title = "Chromosomal aberrations (" + count + ")";
        String content = joiner.length() > 0 ? joiner.toString() : ReportResources.NONE;

        return report()
                .title(
                        cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE),
                        cmp.text(content).setStyle(OrangeFonts.TABLE_CONTENT_STYLE_UNPADDED)
                )
                .setDataSource(new JREmptyDataSource());
    }

    private JasperReportBuilder buildPharmacogenetics()
    {
        String title = "Pharmacogenetics (" + mData.peach.size() + ")";
        return PharmacogeneticsTable.build(title, mData.peach, null);
    }
}
