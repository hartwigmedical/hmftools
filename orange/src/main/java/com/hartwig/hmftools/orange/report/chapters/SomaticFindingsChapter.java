package com.hartwig.hmftools.orange.report.chapters;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.pdfdata.SomaticFindingsData;
import com.hartwig.hmftools.orange.report.tables.ChrArmCopyNumberTable;
import com.hartwig.hmftools.orange.report.tables.DisruptionTable;
import com.hartwig.hmftools.orange.report.tables.DnaFusionTable;
import com.hartwig.hmftools.orange.report.tables.GainDeletionTable;
import com.hartwig.hmftools.orange.report.tables.SignatureAllocationTable;
import com.hartwig.hmftools.orange.report.tables.SomaticVariantTable;
import com.hartwig.hmftools.orange.report.tables.ViralPresenceTable;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Images;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.HorizontalAlignment;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class SomaticFindingsChapter implements ReportChapter
{
    private final SomaticFindingsData mData;
    private final ReportResources mReportResources;

    public SomaticFindingsChapter(final SomaticFindingsData data, final ReportResources reportResources)
    {
        mData = data;
        mReportResources = reportResources;
    }

    @NotNull
    @Override
    public String name()
    {
        return "Somatic Findings";
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

        addSomaticVariants(document);
        addSomaticAmpDels(document);
        addFusions(document);
        addBreakendDisruptions(document);

        if(!mData.tumorOnlyMode)
        {
            addViralPresence(document);
        }

        addChrArmCopyNumbers(document);

        if(!mData.tumorOnlyMode)
        {
            addSignatureAllocations(document);
        }

        addLinxPlots(document);
    }

    private void addSomaticVariants(final Document document)
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

        String titleDrivers = "Small Variants (" + variantImpactCount + ")";

        document.add(SomaticVariantTable.build(
                titleDrivers, contentWidth(), mData.somaticVariants, mReportResources,
                mData.tumorOnlyMode, mData.hasRnaSample));
    }

    private void addSomaticAmpDels(final Document document)
    {
        String titleDrivers = "Amplifications and Deletions (" + mData.somaticGainsDels.size() + ")";

        document.add(GainDeletionTable.build(
                titleDrivers, contentWidth(), mData.somaticGainsDels, mReportResources, mData.hasRna));
    }

    private void addChrArmCopyNumbers(final Document document)
    {
        String title = "Arm Copy Number Abberations";

        document.add(ChrArmCopyNumberTable.build(title, contentWidth(), mData.armCopyNumberAbberations, mReportResources));
    }

    private void addFusions(final Document document)
    {
        String titleDrivers = "Fusions (" + mData.fusions.size() + ")";
        document.add(DnaFusionTable.build(titleDrivers, contentWidth(), mData.fusions, mReportResources));
    }

    private void addViralPresence(final Document document)
    {
        VirusInterpreterData virusInterpreter = mData.virusInterpreter;

        if(virusInterpreter != null)
        {
            String titleDrivers = "Viruses (" + virusInterpreter.reportableViruses().size() + ")";
            document.add(ViralPresenceTable.build(titleDrivers, contentWidth(), virusInterpreter.reportableViruses(), mReportResources));
        }
    }

    private void addBreakendDisruptions(final Document document)
    {
        List<LinxBreakend> somaticBreakends = mData.somaticBreakends;

        Set<Integer> uniqueSvs = somaticBreakends.stream().map(LinxBreakend::svId).collect(Collectors.toSet());

        String titleDriver = "Disruptions (" + uniqueSvs.size() + ")";
        document.add(DisruptionTable.build(titleDriver, contentWidth(), somaticBreakends, mReportResources));
    }

    private void addSignatureAllocations(final Document document)
    {
        if(mData.sigAllocations != null)
        {
            String title = "Signature Allocations (" + mData.sigAllocations.size() + ")";
            document.add(SignatureAllocationTable.build(title, contentWidth(), mData.sigAllocations, mReportResources));
        }
    }

    private void addLinxPlots(final Document document)
    {
        String title = "Structural Driver Plots (" + mData.linxDriverPlotPaths.size() + ")";

        document.add(new Paragraph(title).addStyle(mReportResources.tableTitleStyle()));

        Table table = new Table(2);
        Cells cells = new Cells(mReportResources);

        for(String plotPath : mData.linxDriverPlotPaths)
        {
            Image image = Images.build(plotPath);
            image.setMaxWidth(Math.round(contentWidth() / 2D) - 2);
            image.setHorizontalAlignment(HorizontalAlignment.CENTER);
            table.addCell(cells.createImage(image));
        }

        if(mData.linxDriverPlotPaths.size() % 2 == 1)
        {
            table.addCell(cells.createContent(Strings.EMPTY));
        }

        document.add(table);
    }
}
