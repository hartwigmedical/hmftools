package com.hartwig.hmftools.orange.report.chapters;

import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.algo.QcStatusInterpretation;
import com.hartwig.hmftools.datamodel.sigs.SignatureAllocation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.orange.report.DocumentContext;
import com.hartwig.hmftools.orange.report.PlotPathResolver;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.tables.ChrArmCopyNumberTable;
import com.hartwig.hmftools.orange.report.tables.DisruptionTable;
import com.hartwig.hmftools.orange.report.tables.DnaFusionTable;
import com.hartwig.hmftools.orange.report.tables.GainDeletionTable;
import com.hartwig.hmftools.orange.report.tables.SignatureAllocationTable;
import com.hartwig.hmftools.orange.report.tables.SomaticVariantTable;
import com.hartwig.hmftools.orange.report.tables.ViralPresenceTable;

import org.apache.pdfbox.pdmodel.common.PDRectangle;
import org.jetbrains.annotations.NotNull;

public class SomaticFindingsChapter implements ReportChapter
{
    private final OrangeConfig mConfig;
    private final OrangeRecord mReport;
    private final PlotPathResolver mPlotPathResolver;
    private final ReportResources mReportResources;

    public SomaticFindingsChapter(
            final OrangeConfig config, final OrangeRecord report, final PlotPathResolver plotPathResolver,
            final ReportResources reportResources)
    {
        mConfig = config;
        mReport = report;
        mPlotPathResolver = plotPathResolver;
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

        addSomaticVariants(document);
        document.addSpacing(10);
        addSomaticAmpDels(document);
        document.addSpacing(10);
        addFusions(document);
        document.addSpacing(10);
        addBreakendDisruptions(document);
        document.addSpacing(10);

        if(!mReport.tumorOnlyMode())
        {
            addViralPresence(document);
            document.addSpacing(10);
        }

        addChrArmCopyNumbers(document);
        document.addSpacing(10);

        if(!mReport.tumorOnlyMode())
        {
            addSignatureAllocations(document);
            document.addSpacing(10);
        }

        addLinxPlots(document);
    }

    private void addSomaticVariants(final DocumentContext document) throws IOException
    {
        String driverVariantsTitle = "Small Variants";

        int variantImpactCount = 0;

        for(PurpleVariant variant : mReport.purple().somaticVariants())
        {
            if(variant.canonicalImpact().reported())
            {
                ++variantImpactCount;
            }

            variantImpactCount += variant.otherImpacts().stream().filter(PurpleTranscriptImpact::reported).mapToInt(x -> 1).sum();
        }

        String titleDrivers = driverVariantsTitle + " (" + variantImpactCount + ")";

        document.addTable(SomaticVariantTable.build(
                document, titleDrivers, contentWidth(), mReport.purple().somaticVariants(), mReportResources,
                mReport.tumorOnlyMode(), mConfig != null && mConfig.RnaSampleId != null));
    }

    private void addSomaticAmpDels(final DocumentContext document) throws IOException
    {
        String driverAmpsDelsTitle = "Amplifications and Deletions";
        String titleDrivers = driverAmpsDelsTitle + " (" + mReport.purple().somaticGainsDels().size() + ")";

        document.addTable(GainDeletionTable.build(
                document, titleDrivers, contentWidth(), mReport.purple().somaticGainsDels(), mReportResources, mReport.hasRna()));
    }

    private void addChrArmCopyNumbers(final DocumentContext document) throws IOException
    {
        String title = "Arm Copy Number Aberrations";
        document.addTable(ChrArmCopyNumberTable.build(document, title, contentWidth(), mReport.purple()
                .armCopyNumberAbberations(), mReportResources));
    }

    private void addFusions(final DocumentContext document) throws IOException
    {
        String driverFusionsTitle = "Fusions";
        String titleDrivers = driverFusionsTitle + " (" + mReport.linx().fusions().size() + ")";
        document.addTable(DnaFusionTable.build(document, titleDrivers, contentWidth(), mReport.linx().fusions(), mReportResources));
    }

    private void addViralPresence(final DocumentContext document) throws IOException
    {
        VirusInterpreterData virusInterpreter = mReport.virusInterpreter();

        if(virusInterpreter != null)
        {
            String driverVirusTitle = "Viruses";
            String titleDrivers = driverVirusTitle + " (" + virusInterpreter.reportableViruses().size() + ")";
            document.addTable(ViralPresenceTable.build(document, titleDrivers, contentWidth(), virusInterpreter.reportableViruses(), mReportResources));
        }
    }

    private void addBreakendDisruptions(final DocumentContext document) throws IOException
    {
        String driverGeneDisruptionsTitle = "Disruptions";

        List<LinxBreakend> somaticBreakends = mReport.linx().somaticBreakends();
        Set<Integer> uniqueSvs = somaticBreakends.stream().map(LinxBreakend::svId).collect(Collectors.toSet());

        String titleDriver = driverGeneDisruptionsTitle + " (" + uniqueSvs.size() + ")";
        document.addTable(DisruptionTable.build(document, titleDriver, contentWidth(), somaticBreakends, mReportResources));
    }

    private void addSignatureAllocations(final DocumentContext document) throws IOException
    {
        List<SignatureAllocation> sigAllocations = mReport.sigAllocations();

        if(sigAllocations != null)
        {
            String signatureTitle = "Signature Allocations";
            String title = signatureTitle + " (" + sigAllocations.size() + ")";
            document.addTable(SignatureAllocationTable.build(document, title, contentWidth(), sigAllocations, mReportResources));
        }
    }

    private void addLinxPlots(final DocumentContext document) throws IOException
    {
        String title = "Structural Driver Plots (" + mReport.plots().linxDriverPlots().size() + ")";
        document.addParagraph(title, mReportResources.tableTitleStyle());

        // 2-column grid layout
        float colWidth = contentWidth() / 2f;
        float maxPlotHeight = Math.round(contentWidth() / 2f) - 2;
        float marginLeft = document.marginLeft();

        java.util.List<String> plots = mReport.plots().linxDriverPlots();
        for(int i = 0; i < plots.size(); i += 2)
        {
            float rowY = document.cursorY();

            // Check if we need a new page
            if(rowY - maxPlotHeight < document.contentEndY())
            {
                document.newPage();
                rowY = document.cursorY();
            }

            // Left image
            document.addImageAt(mPlotPathResolver.resolve(plots.get(i)), marginLeft, rowY, colWidth, maxPlotHeight);

            // Right image (if exists)
            if(i + 1 < plots.size())
            {
                document.addImageAt(mPlotPathResolver.resolve(plots.get(i + 1)), marginLeft + colWidth, rowY, colWidth, maxPlotHeight);
            }

            document.setCursorY(rowY - maxPlotHeight - 5);
        }
    }
}
