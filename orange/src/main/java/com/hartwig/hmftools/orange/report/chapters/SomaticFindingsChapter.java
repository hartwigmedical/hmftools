package com.hartwig.hmftools.orange.report.chapters;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.pdfdata.SomaticFindingsData;
import com.hartwig.hmftools.orange.report.tables.ChrArmCopyNumberTable;
import com.hartwig.hmftools.orange.report.tables.DisruptionTable;
import com.hartwig.hmftools.orange.report.tables.DnaFusionTable;
import com.hartwig.hmftools.orange.report.tables.GainDeletionTable;
import com.hartwig.hmftools.orange.report.tables.SignatureAllocationTable;
import com.hartwig.hmftools.orange.report.tables.SomaticVariantTable;
import com.hartwig.hmftools.orange.report.tables.ViralPresenceTable;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;
import net.sf.dynamicreports.report.constant.HorizontalImageAlignment;
import net.sf.dynamicreports.report.constant.PageOrientation;
import net.sf.dynamicreports.report.constant.PageType;
import net.sf.dynamicreports.report.constant.SplitType;
import net.sf.jasperreports.engine.JREmptyDataSource;

public class SomaticFindingsChapter implements ReportChapter
{
    private final SomaticFindingsData mData;

    public SomaticFindingsChapter(final SomaticFindingsData data, final Object unused)
    {
        mData = data;
    }

    @Override
    public String name()
    {
        return "Somatic Findings";
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
        report = report.title(cmp.text(name()).setStyle(OrangeFonts.CHAPTER_TITLE_STYLE));
        VerticalListBuilder content = cmp.verticalList();

        //        content.add(cmp.text(name()).setStyle(OrangeFonts.CHAPTER_TITLE_STYLE));

        if(mData.hasPurpleFail)
        {
            content.add(cmp.text(ReportResources.NOT_AVAILABLE).setStyle(OrangeFonts.TABLE_CONTENT_STYLE));
            return report.summary(content);
        }

        content.add(cmp.subreport(buildSomaticVariants()));
        content.add(cmp.subreport(buildSomaticAmpDels()));
        content.add(cmp.subreport(buildFusions()));
        content.add(cmp.subreport(buildDisruptions()));

        if(!mData.tumorOnlyMode)
        {
            content.add(cmp.subreport(buildViralPresence()));
        }

        content.add(cmp.subreport(buildChrArmCopyNumbers()));

        if(!mData.tumorOnlyMode && mData.sigAllocations != null)
        {
            content.add(cmp.subreport(buildSignatureAllocations()));
        }

        content.add(cmp.pageBreak());
        addLinxPlots(content);

        return report.setSummarySplitType(SplitType.IMMEDIATE).summary(content);
    }

    private JasperReportBuilder buildSomaticVariants()
    {
        int variantImpactCount = 0;
        for(PurpleVariant variant : mData.somaticVariants)
        {
            if(variant.canonicalImpact().reported())
            {
                ++variantImpactCount;
            }
            variantImpactCount += variant.otherImpacts().stream().filter(PurpleTranscriptImpact::reported).mapToInt(x -> 1).sum();
        }
        String title = "Small Variants (" + variantImpactCount + ")";
        return SomaticVariantTable.build(title, mData.somaticVariants, null, mData.tumorOnlyMode, mData.hasRnaSample);
    }

    private JasperReportBuilder buildSomaticAmpDels()
    {
        String title = "Amplifications and Deletions (" + mData.somaticGainsDels.size() + ")";
        return GainDeletionTable.build(title, mData.somaticGainsDels, null, mData.hasRna);
    }

    private JasperReportBuilder buildFusions()
    {
        String title = "Fusions (" + mData.fusions.size() + ")";
        return DnaFusionTable.build(title, mData.fusions, null);
    }

    private JasperReportBuilder buildDisruptions()
    {
        List<LinxBreakend> somaticBreakends = mData.somaticBreakends;
        Set<Integer> uniqueSvs = somaticBreakends.stream().map(LinxBreakend::svId).collect(Collectors.toSet());
        String title = "Disruptions (" + uniqueSvs.size() + ")";
        return DisruptionTable.build(title, somaticBreakends, null);
    }

    private JasperReportBuilder buildViralPresence()
    {
        VirusInterpreterData virusInterpreter = mData.virusInterpreter;
        int count = virusInterpreter != null ? virusInterpreter.reportableViruses().size() : 0;
        String title = "Viruses (" + count + ")";
        List<com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry> viruses =
                virusInterpreter != null ? virusInterpreter.reportableViruses() : java.util.List.of();
        return ViralPresenceTable.build(title, viruses, null);
    }

    private JasperReportBuilder buildChrArmCopyNumbers()
    {
        String title = "Arm Copy Number Abberations";
        return ChrArmCopyNumberTable.build(title, mData.armCopyNumberAbberations, null);
    }

    private JasperReportBuilder buildSignatureAllocations()
    {
        String title = "Signature Allocations (" + mData.sigAllocations.size() + ")";
        return SignatureAllocationTable.build(title, mData.sigAllocations, null);
    }

    private void addLinxPlots(final VerticalListBuilder content)
    {
        if(mData.linxDriverPlotPaths.isEmpty())
        {
            return;
        }

        String title = "Structural Driver Plots (" + mData.linxDriverPlotPaths.size() + ")";
        content.add(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE));

        List<String> paths = mData.linxDriverPlotPaths;
        for(int i = 0; i < paths.size(); i += 2)
        {
            if(i + 1 < paths.size())
            {
                content.add(cmp.horizontalList(
                        cmp.image(paths.get(i))
                                .setFixedHeight(350)
                                .setHorizontalImageAlignment(HorizontalImageAlignment.CENTER),
                        cmp.image(paths.get(i + 1))
                                .setFixedHeight(350)
                                .setHorizontalImageAlignment(HorizontalImageAlignment.CENTER)
                ));
            }
            else
            {
                content.add(cmp.horizontalList(
                        cmp.image(paths.get(i))
                                .setFixedHeight(350)
                                .setHorizontalImageAlignment(HorizontalImageAlignment.CENTER),
                        cmp.text("")
                ));
            }
        }
    }
}
