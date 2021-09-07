package com.hartwig.hmftools.orange.report.chapters;

import java.util.List;

import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.algo.selection.CopyNumberSelector;
import com.hartwig.hmftools.orange.algo.selection.FusionSelector;
import com.hartwig.hmftools.orange.algo.selection.SomaticVariantSelector;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.tables.FusionTable;
import com.hartwig.hmftools.orange.report.tables.GeneCopyNumberTable;
import com.hartwig.hmftools.orange.report.tables.GeneDisruptionTable;
import com.hartwig.hmftools.orange.report.tables.HomozygousDisruptionTable;
import com.hartwig.hmftools.orange.report.tables.SomaticVariantTable;
import com.hartwig.hmftools.orange.report.tables.StructuralDriverTable;
import com.hartwig.hmftools.orange.report.tables.ViralPresenceTable;
import com.hartwig.hmftools.orange.report.util.ImageUtil;
import com.hartwig.hmftools.orange.report.util.TableUtil;
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

    public SomaticFindingsChapter(@NotNull final OrangeReport report) {
        this.report = report;
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
        addGeneDisruptions(document);
        addStructuralDrivers(document);
        addStructuralDriverPlots(document);
    }

    private void addSomaticVariants(@NotNull Document document) {
        String titleDrivers = "Driver variants (" + report.purple().reportableSomaticVariants().size() + ")";
        document.add(SomaticVariantTable.build(titleDrivers, contentWidth(), report.purple().reportableSomaticVariants()));

        List<ReportableVariant> nonDriverVariants = SomaticVariantSelector.selectNonDrivers(report.purple().unreportedSomaticVariants(),
                report.purple().reportableSomaticVariants(),
                report.protect());
        String titleNonDrivers = "Other potentially relevant variants (" + nonDriverVariants.size() + ")";
        document.add(SomaticVariantTable.build(titleNonDrivers, contentWidth(), max10(nonDriverVariants)));
    }

    private void addKataegisPlot(@NotNull Document document) {
        document.add(new Paragraph("Kataegis plot").addStyle(ReportResources.tableTitleStyle()));
        if (report.plots().purpleKataegisPlot() != null) {
            Image image = ImageUtil.build(report.plots().purpleKataegisPlot());
            image.setMaxWidth(contentWidth());
            image.setHorizontalAlignment(HorizontalAlignment.CENTER);
            document.add(image);
        } else {
            document.add(new Paragraph("No kataegis plot could be generated for this sample").addStyle(ReportResources.tableContentStyle()));
        }
    }

    private void addSomaticAmpDels(@NotNull Document document) {
        String titleDrivers = "Driver amps/dels (" + report.purple().reportableGainsLosses().size() + ")";
        document.add(GeneCopyNumberTable.build(titleDrivers, contentWidth(), report.purple().reportableGainsLosses()));

        List<ReportableGainLoss> gains = CopyNumberSelector.selectNonDriverGains(report.purple().unreportedGainsLosses());
        String titleGains = "Other regions with amps (" + gains.size() + ")";
        document.add(GeneCopyNumberTable.build(titleGains, contentWidth(), max10(gains)));

        List<ReportableGainLoss> losses =
                CopyNumberSelector.selectNonDriverLosses(report.purple().unreportedGainsLosses(), report.purple().reportableGainsLosses());
        String titleLosses = "Other autosomal regions with dels (" + losses.size() + ")";
        document.add(GeneCopyNumberTable.build(titleLosses, contentWidth(), max10(losses)));
    }

    private void addFusions(@NotNull Document document) {
        String titleDrivers = "Driver fusions (" + report.linx().reportableFusions().size() + ")";
        document.add(FusionTable.build(titleDrivers, contentWidth(), report.linx().reportableFusions()));

        List<LinxFusion> nonDriverFusions = FusionSelector.selectNonDriverFusions(report.linx().unreportedFusions(), report.protect());
        String titleNonDrivers = "Other potentially interesting fusions (" + nonDriverFusions.size() + ")";
        document.add(FusionTable.build(titleNonDrivers, contentWidth(), max10(nonDriverFusions)));
    }

    private void addViralPresence(@NotNull Document document) {
        String titleDrivers = "Driver viruses (" + report.virusInterpreter().reportableViruses().size() + ")";
        document.add(ViralPresenceTable.build(titleDrivers, contentWidth(), report.virusInterpreter().reportableViruses()));

        String titleNonDrivers = "Other viral presence (" + report.virusInterpreter().unreportedViruses().size() + ")";
        document.add(ViralPresenceTable.build(titleNonDrivers, contentWidth(), report.virusInterpreter().unreportedViruses()));
    }

    private void addHomozygousDisruptions(@NotNull Document document) {
        String titleDrivers = "Homozygous disruptions (" + report.linx().homozygousDisruptions().size() + ")";
        document.add(HomozygousDisruptionTable.build(titleDrivers, contentWidth(), report.linx().homozygousDisruptions()));
    }

    private void addGeneDisruptions(@NotNull Document document) {
        String titleDrivers = "Gene disruptions (" + report.linx().geneDisruptions().size() + ")";
        document.add(GeneDisruptionTable.build(titleDrivers, contentWidth(), report.linx().geneDisruptions()));
    }

    private void addStructuralDrivers(final Document document) {
        String titleDrivers = "Structural drivers (" + report.linx().drivers().size() + ")";
        document.add(StructuralDriverTable.build(titleDrivers, contentWidth(), report.linx().drivers()));
    }

    private void addStructuralDriverPlots(@NotNull Document document) {
        String title = "Structural driver plots (" + report.plots().linxDriverPlots().size() + ")";
        document.add(new Paragraph(title).addStyle(ReportResources.tableTitleStyle()));
        Table table = new Table(2);
        for (String plot : report.plots().linxDriverPlots()) {
            Image image = ImageUtil.build(plot);
            image.setMaxWidth(Math.round(contentWidth() / 2D) - 2);
            image.setHorizontalAlignment(HorizontalAlignment.CENTER);
            table.addCell(TableUtil.createImageCell(image));
        }

        if (report.plots().linxDriverPlots().size() % 2 == 1) {
            table.addCell(TableUtil.createContentCell(Strings.EMPTY));
        }

        document.add(table);
    }

    @NotNull
    private static <T> List<T> max10(@NotNull List<T> elements) {
        return elements.subList(0, Math.min(10, elements.size()));
    }
}
