package com.hartwig.hmftools.orange.report.chapters;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleGainLoss;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.sigs.SignatureAllocation;
import com.hartwig.hmftools.datamodel.virus.AnnotatedVirus;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.orange.report.PlotPathResolver;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.datamodel.BreakendEntry;
import com.hartwig.hmftools.orange.report.datamodel.BreakendEntryFactory;
import com.hartwig.hmftools.orange.report.datamodel.VariantEntry;
import com.hartwig.hmftools.orange.report.datamodel.VariantEntryFactory;
import com.hartwig.hmftools.orange.report.interpretation.VariantDedup;
import com.hartwig.hmftools.orange.report.tables.BreakendTable;
import com.hartwig.hmftools.orange.report.tables.DNAFusionTable;
import com.hartwig.hmftools.orange.report.tables.GainLossTable;
import com.hartwig.hmftools.orange.report.tables.HomozygousDisruptionTable;
import com.hartwig.hmftools.orange.report.tables.LossOfHeterozygosityTable;
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

public class SomaticFindingsChapter implements ReportChapter {

    @NotNull
    private final OrangeRecord report;
    @NotNull
    private final PlotPathResolver plotPathResolver;

    public SomaticFindingsChapter(@NotNull final OrangeRecord report, @NotNull final PlotPathResolver plotPathResolver) {
        this.report = report;
        this.plotPathResolver = plotPathResolver;
    }

    @NotNull
    @Override
    public String name() {
        return "Somatic Findings";
    }

    @NotNull
    @Override
    public PageSize pageSize() {
        return PageSize.A4;
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph(name()).addStyle(ReportResources.chapterTitleStyle()));

        addSomaticVariants(document);
        addKataegisPlot(document);
        addSomaticAmpDels(document);
        addFusions(document);
        addViralPresence(document);
        addHomozygousDisruptions(document);
        addBreakends(document);
        addLossOfHeterozygosity(document);
        addSignatureAllocations(document);
        addStructuralDriverPlots(document);
    }

    private void addSomaticVariants(@NotNull Document document) {
        List<PurpleDriver> somaticDrivers = report.purple().somaticDrivers();

        List<VariantEntry> reportableVariants =
                VariantEntryFactory.create(VariantDedup.apply(report.purple().reportableSomaticVariants()), somaticDrivers);
        String titleDrivers = "Driver variants (" + reportableVariants.size() + ")";
        document.add(SomaticVariantTable.build(titleDrivers, contentWidth(), reportableVariants));

        List<VariantEntry> additionalSuspectVariants =
                VariantEntryFactory.create(VariantDedup.apply(report.purple().additionalSuspectSomaticVariants()), somaticDrivers);
        String titleNonDrivers = "Other potentially relevant variants (" + additionalSuspectVariants.size() + ")";
        document.add(SomaticVariantTable.build(titleNonDrivers, contentWidth(), max10(additionalSuspectVariants)));
    }

    private void addKataegisPlot(@NotNull Document document) {
        document.add(new Paragraph("Kataegis plot").addStyle(ReportResources.tableTitleStyle()));
        String kataegisPlot = report.plots().purpleKataegisPlot();
        if (kataegisPlot != null) {
            Image image = Images.build(plotPathResolver.resolve(kataegisPlot));
            image.setMaxWidth(contentWidth());
            image.setHorizontalAlignment(HorizontalAlignment.CENTER);
            document.add(image);
        } else {
            document.add(new Paragraph("No kataegis plot could be generated for this sample").addStyle(ReportResources.tableContentStyle()));
        }
    }

    private void addSomaticAmpDels(@NotNull Document document) {
        String titleDrivers = "Driver amps/dels (" + report.purple().reportableSomaticGainsLosses().size() + ")";
        document.add(GainLossTable.build(titleDrivers,
                contentWidth(),
                report.purple().reportableSomaticGainsLosses(),
                report.isofox()));

        String titleNearDriverGains =
                "Potentially interesting near-driver amps (" + report.purple().nearReportableSomaticGains().size() + ")";
        document.add(GainLossTable.build(titleNearDriverGains,
                contentWidth(),
                report.purple().nearReportableSomaticGains(),
                report.isofox()));

        List<PurpleGainLoss> suspectGains = selectGains(report.purple().additionalSuspectSomaticGainsLosses());
        String titleSuspectGains = "Other regions with amps (" + suspectGains.size() + ")";
        document.add(GainLossTable.build(titleSuspectGains, contentWidth(), max10(suspectGains), report.isofox()));

        List<PurpleGainLoss> suspectLosses = selectLosses(report.purple().additionalSuspectSomaticGainsLosses());
        String titleSuspectLosses = "Regions with deletions in genes in other autosomal regions (" + suspectLosses.size() + ")";
        document.add(GainLossTable.build(titleSuspectLosses, contentWidth(), max10(suspectLosses), report.isofox()));
    }

    @NotNull
    private static List<PurpleGainLoss> selectGains(@NotNull List<PurpleGainLoss> gainsLosses) {
        List<PurpleGainLoss> gains = Lists.newArrayList();
        for (PurpleGainLoss gainLoss : gainsLosses) {
            if (gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_GAIN
                    || gainLoss.interpretation() == CopyNumberInterpretation.FULL_GAIN) {
                gains.add(gainLoss);
            }
        }
        return gains;
    }

    @NotNull
    private static List<PurpleGainLoss> selectLosses(@NotNull List<PurpleGainLoss> gainsLosses) {
        List<PurpleGainLoss> losses = Lists.newArrayList();
        for (PurpleGainLoss gainLoss : gainsLosses) {
            if (gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_LOSS
                    || gainLoss.interpretation() == CopyNumberInterpretation.FULL_LOSS) {
                losses.add(gainLoss);
            }
        }
        return losses;
    }

    private void addFusions(@NotNull Document document) {
        String titleDrivers = "Driver fusions (" + report.linx().reportableSomaticFusions().size() + ")";
        document.add(DNAFusionTable.build(titleDrivers, contentWidth(), report.linx().reportableSomaticFusions(), report.isofox()));

        String titleNonDrivers = "Other potentially interesting fusions (" + report.linx().additionalSuspectSomaticFusions().size() + ")";
        document.add(DNAFusionTable.build(titleNonDrivers,
                contentWidth(),
                max10(report.linx().additionalSuspectSomaticFusions()),
                report.isofox()));
    }

    private void addViralPresence(@NotNull Document document) {
        VirusInterpreterData virusInterpreter = report.virusInterpreter();

        if (virusInterpreter != null) {
            String titleDrivers = "Driver viruses (" + virusInterpreter.reportableViruses().size() + ")";
            document.add(ViralPresenceTable.build(titleDrivers, contentWidth(), virusInterpreter.reportableViruses()));

            List<AnnotatedVirus> unreported = Lists.newArrayList();
            for (AnnotatedVirus virus : virusInterpreter.allViruses()) {
                if (!virus.reported()) {
                    unreported.add(virus);
                }
            }

            String titleNonDrivers = "Other viral presence (" + unreported.size() + ")";
            document.add(ViralPresenceTable.build(titleNonDrivers, contentWidth(), unreported));
        }
    }

    private void addHomozygousDisruptions(@NotNull Document document) {
        String title = "Homozygous disruptions (" + report.linx().somaticHomozygousDisruptions().size() + ")";
        document.add(HomozygousDisruptionTable.build(title, contentWidth(), report.linx().somaticHomozygousDisruptions()));
    }

    private void addBreakends(@NotNull Document document) {
        List<BreakendEntry> reportableBreakends =
                BreakendEntryFactory.create(report.linx().reportableSomaticBreakends(), report.linx().allSomaticStructuralVariants());

        String titleDriver = "Driver gene disruptions (" + reportableBreakends.size() + ")";
        document.add(BreakendTable.build(titleDriver, contentWidth(), reportableBreakends));

        List<BreakendEntry> additionalSuspectBreakends =
                BreakendEntryFactory.create(report.linx().additionalSuspectSomaticBreakends(), report.linx().allSomaticStructuralVariants());
        String titleNonDrivers = "Other potentially interesting gene disruptions (" + additionalSuspectBreakends.size() + ")";
        document.add(BreakendTable.build(titleNonDrivers, contentWidth(), additionalSuspectBreakends));
    }

    private void addLossOfHeterozygosity(@NotNull Document document) {
        List<PurpleGeneCopyNumber> suspectGeneCopyNumbersWithLOH = report.purple().suspectGeneCopyNumbersWithLOH();
        String title = "Potentially interesting LOH events in case of MSI or HRD (" + suspectGeneCopyNumbersWithLOH.size() + ")";
        document.add(LossOfHeterozygosityTable.build(title, contentWidth(), suspectGeneCopyNumbersWithLOH));
    }

    private void addSignatureAllocations(@NotNull Document document) {
        List<SignatureAllocation> sigAllocations = report.sigAllocations();

        if (sigAllocations != null) {
            String title = "Signature allocations (" + sigAllocations.size() + ")";
            document.add(SignatureAllocationTable.build(title, contentWidth(), sigAllocations));
        }
    }

    private void addStructuralDriverPlots(@NotNull Document document) {
        String title = "Structural driver plots (" + report.plots().linxDriverPlots().size() + ")";
        document.add(new Paragraph(title).addStyle(ReportResources.tableTitleStyle()));
        Table table = new Table(2);
        for (String plot : report.plots().linxDriverPlots()) {
            Image image = Images.build(plotPathResolver.resolve(plot));
            image.setMaxWidth(Math.round(contentWidth() / 2D) - 2);
            image.setHorizontalAlignment(HorizontalAlignment.CENTER);
            table.addCell(Cells.createImage(image));
        }

        if (report.plots().linxDriverPlots().size() % 2 == 1) {
            table.addCell(Cells.createContent(Strings.EMPTY));
        }

        document.add(table);
    }

    @NotNull
    private static <T> List<T> max10(@NotNull List<T> elements) {
        return elements.subList(0, Math.min(10, elements.size()));
    }
}
