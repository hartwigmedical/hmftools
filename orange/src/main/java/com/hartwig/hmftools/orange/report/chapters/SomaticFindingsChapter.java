package com.hartwig.hmftools.orange.report.chapters;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.linx.FusionLikelihoodType;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleGainLoss;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.sigs.SignatureAllocation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;
import com.hartwig.hmftools.orange.algo.purple.DriverInterpretation;
import com.hartwig.hmftools.orange.report.PlotPathResolver;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.datamodel.BreakendEntry;
import com.hartwig.hmftools.orange.report.datamodel.BreakendEntryFactory;
import com.hartwig.hmftools.orange.report.datamodel.VariantEntry;
import com.hartwig.hmftools.orange.report.datamodel.VariantEntryFactory;
import com.hartwig.hmftools.orange.report.interpretation.PurpleQCInterpretation;
import com.hartwig.hmftools.orange.report.interpretation.VariantDedup;
import com.hartwig.hmftools.orange.report.tables.BreakendTable;
import com.hartwig.hmftools.orange.report.tables.DnaFusionTable;
import com.hartwig.hmftools.orange.report.tables.GainLossTable;
import com.hartwig.hmftools.orange.report.tables.HomozygousDisruptionTable;
import com.hartwig.hmftools.orange.report.tables.LossOfHeterozygosityTable;
import com.hartwig.hmftools.orange.report.tables.SignatureAllocationTable;
import com.hartwig.hmftools.orange.report.tables.SomaticVariantTable;
import com.hartwig.hmftools.orange.report.tables.ViralPresenceTable;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Images;
import com.hartwig.hmftools.orange.report.util.Tables;
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
    @NotNull
    private final OrangeRecord report;
    @NotNull
    private final PlotPathResolver plotPathResolver;
    @NotNull
    private final ReportResources reportResources;

    public SomaticFindingsChapter(@NotNull final OrangeRecord report, @NotNull final PlotPathResolver plotPathResolver,
            @NotNull final ReportResources reportResources)
    {
        this.report = report;
        this.plotPathResolver = plotPathResolver;
        this.reportResources = reportResources;
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
        document.add(new Paragraph(name()).addStyle(reportResources.chapterTitleStyle()));

        addSomaticVariants(document);
        if(!PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            addKataegisPlot(document);
        }
        addSomaticAmpDels(document);
        addFusions(document);

        if(!report.tumorOnlyMode())
        {
            addViralPresence(document);
        }

        addHomozygousDisruptions(document);
        addBreakends(document);
        addLossOfHeterozygosity(document);

        if(!report.tumorOnlyMode())
        {
            addSignatureAllocations(document);
        }

        if(!PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            addStructuralDriverPlots(document);
        }
    }

    private void addSomaticVariants(@NotNull Document document)
    {
        String driverVariantsTitle = "Driver variants";
        String otherPotentiallyInterestingTitle = "Other potentially relevant variants";

        if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            Tables tables = new Tables(reportResources);
            document.add(tables.createNotAvailable(driverVariantsTitle, contentWidth()));
            document.add(tables.createNotAvailable(otherPotentiallyInterestingTitle, contentWidth()));
        }
        else
        {
            List<PurpleDriver> somaticDrivers = report.purple().somaticDrivers();

            List<VariantEntry> reportableVariants =
                    VariantEntryFactory.create(VariantDedup.apply(report.purple().reportableSomaticVariants()), somaticDrivers);
            String titleDrivers = driverVariantsTitle + " (" + reportableVariants.size() + ")";
            document.add(SomaticVariantTable.build(titleDrivers, contentWidth(), reportableVariants, reportResources));

            List<VariantEntry> additionalSuspectVariants =
                    VariantEntryFactory.create(VariantDedup.apply(report.purple().additionalSuspectSomaticVariants()), somaticDrivers);
            String titleNonDrivers = otherPotentiallyInterestingTitle + " (" + additionalSuspectVariants.size() + ")";
            document.add(SomaticVariantTable.build(titleNonDrivers, contentWidth(), max10(additionalSuspectVariants), reportResources));
        }
    }

    private void addKataegisPlot(@NotNull Document document)
    {
        document.add(new Paragraph("Kataegis plot").addStyle(reportResources.tableTitleStyle()));
        String kataegisPlot = report.plots().purpleKataegisPlot();
        if(kataegisPlot != null)
        {
            Image image = Images.build(plotPathResolver.resolve(kataegisPlot));
            image.setMaxWidth(contentWidth());
            image.setHorizontalAlignment(HorizontalAlignment.CENTER);
            document.add(image);
        }
        else
        {
            document.add(new Paragraph("No kataegis plot could be generated for this sample")
                    .addStyle(reportResources.tableContentStyle()));
        }
    }

    private void addSomaticAmpDels(@NotNull Document document)
    {
        String driverAmpsDelsTitle = "Driver amps/dels";
        String nearDriverGainsTitle = "Potentially interesting near-driver amps";
        String suspectGainsTitle = "Other regions with amps";
        String suspectLossesTitle = "Regions with deletions in genes in other autosomal regions";

        if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            Tables tables = new Tables(reportResources);
            document.add(tables.createNotAvailable(driverAmpsDelsTitle, contentWidth()));
            document.add(tables.createNotAvailable(nearDriverGainsTitle, contentWidth()));
            document.add(tables.createNotAvailable(suspectGainsTitle, contentWidth()));
            document.add(tables.createNotAvailable(suspectLossesTitle, contentWidth()));
        }
        else
        {
            String titleDrivers = driverAmpsDelsTitle + " (" + report.purple().reportableSomaticGainsLosses().size() + ")";
            document.add(GainLossTable.build(titleDrivers,
                    contentWidth(),
                    report.purple().reportableSomaticGainsLosses(),
                    report.isofox(),
                    reportResources));

            String titleNearDriverGains = nearDriverGainsTitle + " (" + report.purple().nearReportableSomaticGains().size() + ")";
            document.add(GainLossTable.build(titleNearDriverGains,
                    contentWidth(),
                    report.purple().nearReportableSomaticGains(),
                    report.isofox(),
                    reportResources));

            List<PurpleGainLoss> suspectGains = selectGains(report.purple().additionalSuspectSomaticGainsLosses());
            String titleSuspectGains = suspectGainsTitle + " (" + suspectGains.size() + ")";
            document.add(GainLossTable.build(titleSuspectGains, contentWidth(), max10(suspectGains), report.isofox(), reportResources));

            List<PurpleGainLoss> suspectLosses = selectLosses(report.purple().additionalSuspectSomaticGainsLosses());
            String titleSuspectLosses = suspectLossesTitle + " (" + suspectLosses.size() + ")";
            document.add(GainLossTable.build(titleSuspectLosses, contentWidth(), max10(suspectLosses), report.isofox(), reportResources));
        }
    }

    @NotNull
    private static List<PurpleGainLoss> selectGains(@NotNull List<PurpleGainLoss> gainsLosses)
    {
        List<PurpleGainLoss> gains = Lists.newArrayList();
        for(PurpleGainLoss gainLoss : gainsLosses)
        {
            if(gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_GAIN
                    || gainLoss.interpretation() == CopyNumberInterpretation.FULL_GAIN)
            {
                gains.add(gainLoss);
            }
        }
        return gains;
    }

    @NotNull
    private static List<PurpleGainLoss> selectLosses(@NotNull List<PurpleGainLoss> gainsLosses)
    {
        List<PurpleGainLoss> losses = Lists.newArrayList();
        for(PurpleGainLoss gainLoss : gainsLosses)
        {
            if(gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_LOSS
                    || gainLoss.interpretation() == CopyNumberInterpretation.FULL_LOSS)
            {
                losses.add(gainLoss);
            }
        }
        return losses;
    }

    private void addFusions(@NotNull Document document)
    {
        String driverFusionsTitle = "Driver fusions";
        String otherFusionsTitle = "Other potentially interesting fusions";
        String otherFusionsInCaseNoHighDriversTitle = "Potentially interesting in-frame fusions in case no high drivers detected";

        if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            Tables tables = new Tables(reportResources);
            document.add(tables.createNotAvailable(driverFusionsTitle, contentWidth()));
            document.add(tables.createNotAvailable(otherFusionsTitle, contentWidth()));
            document.add(tables.createNotAvailable(otherFusionsInCaseNoHighDriversTitle, contentWidth()));
        }
        else
        {
            String titleDrivers = driverFusionsTitle + " (" + report.linx().reportableSomaticFusions().size() + ")";
            document.add(DnaFusionTable.build(titleDrivers,
                    contentWidth(),
                    report.linx().reportableSomaticFusions(),
                    report.isofox(),
                    reportResources));

            String titleOtherFusions = otherFusionsTitle + " (" + report.linx().additionalSuspectSomaticFusions().size() + ")";
            document.add(DnaFusionTable.build(titleOtherFusions,
                    contentWidth(),
                    max10(report.linx().additionalSuspectSomaticFusions()),
                    report.isofox(),
                    reportResources));

            if(!hasHighDriverEvents(report.linx().allSomaticFusions(), report.purple().somaticDrivers()))
            {
                String titleOtherFusionsNoHighDrivers =
                        otherFusionsInCaseNoHighDriversTitle + " (" + report.linx().additionalViableSomaticFusions().size() + ")";
                document.add(DnaFusionTable.build(titleOtherFusionsNoHighDrivers,
                        contentWidth(),
                        report.linx().additionalViableSomaticFusions(),
                        report.isofox(),
                        reportResources));
            }
            else
            {
                document.add(new Tables(reportResources).createNonContent(otherFusionsInCaseNoHighDriversTitle, contentWidth(),
                        "High driver likelihood events are detected in this sample, therefore this section is empty"));
            }
        }
    }

    private void addViralPresence(@NotNull Document document)
    {
        VirusInterpreterData virusInterpreter = report.virusInterpreter();

        if(virusInterpreter != null)
        {
            String driverVirusTitle = "Driver viruses";
            String nonDriverVirusTitle = "Other viral presence";

            if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
            {
                Tables tables = new Tables(reportResources);
                document.add(tables.createNotAvailable(driverVirusTitle, contentWidth()));
                document.add(tables.createNotAvailable(nonDriverVirusTitle, contentWidth()));
            }
            else
            {
                String titleDrivers = driverVirusTitle + " (" + virusInterpreter.reportableViruses().size() + ")";
                document.add(ViralPresenceTable.build(titleDrivers, contentWidth(), virusInterpreter.reportableViruses(), reportResources));

                List<VirusInterpreterEntry> unreported = Lists.newArrayList();
                for(VirusInterpreterEntry virus : virusInterpreter.allViruses())
                {
                    if(!virus.reported())
                    {
                        unreported.add(virus);
                    }
                }

                String titleNonDrivers = nonDriverVirusTitle + " (" + unreported.size() + ")";
                document.add(ViralPresenceTable.build(titleNonDrivers, contentWidth(), unreported, reportResources));
            }
        }
    }

    private void addHomozygousDisruptions(@NotNull Document document)
    {
        String homozygousDisruptionTitle = "Homozygous disruptions";

        if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            Tables tables = new Tables(reportResources);
            document.add(tables.createNotAvailable(homozygousDisruptionTitle, contentWidth()));
        }
        else
        {
            String title = homozygousDisruptionTitle + " (" + report.linx().somaticHomozygousDisruptions().size() + ")";
            document.add(HomozygousDisruptionTable.build(title,
                    contentWidth(),
                    report.linx().somaticHomozygousDisruptions(),
                    reportResources));
        }
    }

    private void addBreakends(@NotNull Document document)
    {
        String driverGeneDisruptionsTitle = "Driver gene disruptions";
        String nonDriverGeneDisruptionsTitle = "Other potentially interesting gene disruptions";

        if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            Tables tables = new Tables(reportResources);
            document.add(tables.createNotAvailable(driverGeneDisruptionsTitle, contentWidth()));
            document.add(tables.createNotAvailable(nonDriverGeneDisruptionsTitle, contentWidth()));
        }
        else
        {
            List<BreakendEntry> reportableBreakends = BreakendEntryFactory.create(report.linx().reportableSomaticBreakends(),
                    report.linx().allSomaticStructuralVariants(),
                    report.linx().somaticDrivers());

            String titleDriver = driverGeneDisruptionsTitle + " (" + reportableBreakends.size() + ")";
            document.add(BreakendTable.build(titleDriver, contentWidth(), reportableBreakends, reportResources));

            List<BreakendEntry> additionalSuspectBreakends = BreakendEntryFactory.create(report.linx().additionalSuspectSomaticBreakends(),
                    report.linx().allSomaticStructuralVariants(),
                    report.linx().somaticDrivers());
            String titleNonDrivers = nonDriverGeneDisruptionsTitle + " (" + additionalSuspectBreakends.size() + ")";
            document.add(BreakendTable.build(titleNonDrivers, contentWidth(), additionalSuspectBreakends, reportResources));
        }
    }

    private void addLossOfHeterozygosity(@NotNull Document document)
    {
        String lohTitle = "Potentially interesting LOH events";

        if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            Tables tables = new Tables(reportResources);
            document.add(tables.createNotAvailable(lohTitle, contentWidth()));
        }
        else
        {
            List<PurpleGeneCopyNumber> suspectGeneCopyNumbersWithLOH = report.purple().suspectGeneCopyNumbersWithLOH();
            String title = lohTitle + " (" + suspectGeneCopyNumbersWithLOH.size() + ")";
            document.add(LossOfHeterozygosityTable.build(title, contentWidth(), suspectGeneCopyNumbersWithLOH, reportResources));
        }
    }

    private void addSignatureAllocations(@NotNull Document document)
    {
        List<SignatureAllocation> sigAllocations = report.sigAllocations();

        if(sigAllocations != null)
        {
            String signatureTitle = "Signature allocations";

            if(PurpleQCInterpretation.isFail(report.purple().fit().qc()))
            {
                Tables tables = new Tables(reportResources);
                document.add(tables.createNotAvailable(signatureTitle, contentWidth()));
            }
            else
            {
                String title = signatureTitle + " (" + sigAllocations.size() + ")";
                document.add(SignatureAllocationTable.build(title, contentWidth(), sigAllocations, reportResources));
            }
        }
    }

    private void addStructuralDriverPlots(@NotNull Document document)
    {
        String title = "Structural driver plots (" + report.plots().linxDriverPlots().size() + ")";
        document.add(new Paragraph(title).addStyle(reportResources.tableTitleStyle()));
        Table table = new Table(2);
        Cells cells = new Cells(reportResources);
        for(String plot : report.plots().linxDriverPlots())
        {
            Image image = Images.build(plotPathResolver.resolve(plot));
            image.setMaxWidth(Math.round(contentWidth() / 2D) - 2);
            image.setHorizontalAlignment(HorizontalAlignment.CENTER);
            table.addCell(cells.createImage(image));
        }

        if(report.plots().linxDriverPlots().size() % 2 == 1)
        {
            table.addCell(cells.createContent(Strings.EMPTY));
        }

        document.add(table);
    }

    @NotNull
    private static <T> List<T> max10(@NotNull List<T> elements)
    {
        return elements.subList(0, Math.min(10, elements.size()));
    }

    private static boolean hasHighDriverEvents(@NotNull List<LinxFusion> somaticFusions, @NotNull List<PurpleDriver> drivers)
    {
        for(LinxFusion fusion : somaticFusions)
        {
            if(fusion.driverLikelihood() == FusionLikelihoodType.HIGH)
            {
                return true;
            }
        }

        for(PurpleDriver driver : drivers)
        {
            if(DriverInterpretation.interpret(driver.driverLikelihood()) == DriverInterpretation.HIGH)
            {
                return true;
            }
        }

        return false;
    }
}
