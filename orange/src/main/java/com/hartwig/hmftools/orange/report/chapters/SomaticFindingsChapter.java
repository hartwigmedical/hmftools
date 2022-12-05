package com.hartwig.hmftools.orange.report.chapters;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.algo.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.orange.algo.purple.PurpleGainLoss;
import com.hartwig.hmftools.orange.report.PlotPathResolver;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.VariantDedup;
import com.hartwig.hmftools.orange.report.interpretation.VariantEntry;
import com.hartwig.hmftools.orange.report.interpretation.VariantEntryFactory;
import com.hartwig.hmftools.orange.report.tables.BreakendTable;
import com.hartwig.hmftools.orange.report.tables.DNAFusionTable;
import com.hartwig.hmftools.orange.report.tables.GeneCopyNumberTable;
import com.hartwig.hmftools.orange.report.tables.HomozygousDisruptionTable;
import com.hartwig.hmftools.orange.report.tables.LossOfHeterozygosityTable;
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
    private final OrangeReport report;
    @NotNull
    private final PlotPathResolver plotPathResolver;

    public SomaticFindingsChapter(@NotNull final OrangeReport report, @NotNull final PlotPathResolver plotPathResolver) {
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
        addStructuralDriverPlots(document);
    }

    private void addSomaticVariants(@NotNull Document document) {
        List<DriverCatalog> somaticDrivers = report.purple().somaticDrivers();

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
        document.add(GeneCopyNumberTable.build(titleDrivers,
                contentWidth(),
                report.purple().reportableSomaticGainsLosses(),
                report.isofox()));

        String titleNearDriverGains =
                "Potentially interesting near-driver amps (" + report.purple().nearReportableSomaticGains().size() + ")";
        document.add(GeneCopyNumberTable.build(titleNearDriverGains,
                contentWidth(),
                report.purple().nearReportableSomaticGains(),
                report.isofox()));

        List<PurpleGainLoss> suspectGains = selectGains(report.purple().additionalSuspectSomaticGainsLosses());
        String titleSuspectGains = "Other regions with amps (" + suspectGains.size() + ")";
        document.add(GeneCopyNumberTable.build(titleSuspectGains, contentWidth(), max10(suspectGains), report.isofox()));

        List<PurpleGainLoss> suspectLosses = selectLosses(report.purple().additionalSuspectSomaticGainsLosses());
        String titleSuspectLosses = "Regions with deletions in genes in other autosomal regions (" + suspectLosses.size() + ")";
        document.add(GeneCopyNumberTable.build(titleSuspectLosses, contentWidth(), max10(suspectLosses), report.isofox()));
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
        String titleDrivers = "Driver fusions (" + report.linx().reportableFusions().size() + ")";
        document.add(DNAFusionTable.build(titleDrivers, contentWidth(), report.linx().reportableFusions(), report.isofox()));

        String titleNonDrivers = "Other potentially interesting fusions (" + report.linx().additionalSuspectFusions().size() + ")";
        document.add(DNAFusionTable.build(titleNonDrivers,
                contentWidth(),
                max10(report.linx().additionalSuspectFusions()),
                report.isofox()));
    }

    private void addViralPresence(@NotNull Document document) {
        String titleDrivers = "Driver viruses (" + report.virusInterpreter().reportableViruses().size() + ")";
        document.add(ViralPresenceTable.build(titleDrivers, contentWidth(), report.virusInterpreter().reportableViruses()));

        String titleNonDrivers = "Other viral presence (" + report.virusInterpreter().unreportedViruses().size() + ")";
        document.add(ViralPresenceTable.build(titleNonDrivers, contentWidth(), report.virusInterpreter().unreportedViruses()));
    }

    private void addHomozygousDisruptions(@NotNull Document document) {
        String title = "Homozygous disruptions (" + report.linx().homozygousDisruptions().size() + ")";
        document.add(HomozygousDisruptionTable.build(title, contentWidth(), report.linx().homozygousDisruptions()));
    }

    private void addBreakends(@NotNull Document document) {
        List<LinxBreakend> reportableBreakends = report.linx().reportableBreakends();
        String titleDriver = "Driver gene disruptions (" + reportableBreakends.size() + ")";
        document.add(BreakendTable.build(titleDriver, contentWidth(), reportableBreakends));

        List<LinxBreakend> additionalSuspectBreakends = report.linx().additionalSuspectBreakends();
        String titleNonDrivers = "Other potentially interesting gene disruptions (" + additionalSuspectBreakends.size() + ")";
        document.add(BreakendTable.build(titleNonDrivers, contentWidth(), additionalSuspectBreakends));
    }

    private void addLossOfHeterozygosity(@NotNull Document document) {
        List<GeneCopyNumber> suspectGeneCopyNumbersWithLOH = report.purple().suspectGeneCopyNumbersWithLOH();
        String title = "Potentially interesting LOH events in case of MSI or HRD (" + suspectGeneCopyNumbersWithLOH.size() + ")";
        document.add(LossOfHeterozygosityTable.build(title, contentWidth(), suspectGeneCopyNumbersWithLOH));
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
