package com.hartwig.hmftools.orange.report.chapters;

import java.net.MalformedURLException;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.common.linx.ReportableHomozygousDisruption;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.variant.DriverInterpretation;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.components.TableUtil;
import com.itextpdf.io.IOException;
import com.itextpdf.io.image.ImageDataFactory;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.HorizontalAlignment;
import com.itextpdf.layout.property.UnitValue;

import org.jetbrains.annotations.NotNull;

public class FrontPageChapter implements ReportChapter {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#.#");

    @NotNull
    private final OrangeReport report;

    public FrontPageChapter(@NotNull final OrangeReport report) {
        this.report = report;
    }

    @NotNull
    @Override
    public String name() {
        return "Front Page";
    }

    @Override
    public void render(@NotNull Document document) {
        addSummaryTable(document);
        addDetailsAndPlots(document);
    }

    private void addSummaryTable(@NotNull Document document) {
        Table primaryTumorTable = TableUtil.createReportContentTable(new float[] { 1, 1, 1 },
                new Cell[] { TableUtil.createHeaderCell("Configured Primary Tumor"), TableUtil.createHeaderCell("Cuppa Primary Tumor"),
                        TableUtil.createHeaderCell("QC") });
        primaryTumorTable.addCell(TableUtil.createContentCell(toConfiguredTumorType(report.configuredPrimaryTumor())));
        primaryTumorTable.addCell(TableUtil.createContentCell(report.cuppaPrimaryTumor()));
        primaryTumorTable.addCell(TableUtil.createContentCell(purpleQCString()));
        document.add(TableUtil.createWrappingReportTable(primaryTumorTable));
    }

    private void addDetailsAndPlots(@NotNull Document document) {
        Table topTable = new Table(UnitValue.createPercentArray(new float[] { 1, 1 })).setWidth(ReportResources.CONTENT_WIDTH_WIDE - 5);

        Table summary = new Table(UnitValue.createPercentArray(new float[] { 1, 1 }));
        summary.addCell(TableUtil.createKeyCell("Purity:"));
        summary.addCell(TableUtil.createValueCell(purityString()));
        summary.addCell(TableUtil.createKeyCell("Ploidy:"));
        summary.addCell(TableUtil.createValueCell(ploidyString()));
        summary.addCell(TableUtil.createKeyCell("Fit method:"));
        summary.addCell(TableUtil.createValueCell(report.purple().fittedPurityMethod().toString()));
        summary.addCell(TableUtil.createKeyCell("Whole genome duplicated:"));
        summary.addCell(TableUtil.createValueCell(report.purple().wholeGenomeDuplication() ? "Yes" : "No"));
        summary.addCell(TableUtil.createKeyCell("Somatic variant drivers:"));
        summary.addCell(TableUtil.createValueCell(somaticDriverString()));
        summary.addCell(TableUtil.createKeyCell("Germline variant drivers:"));
        summary.addCell(TableUtil.createValueCell(germlineDriverString()));
        summary.addCell(TableUtil.createKeyCell("Copy number drivers:"));
        summary.addCell(TableUtil.createValueCell(copyNumberDriverString()));
        summary.addCell(TableUtil.createKeyCell("Disruption drivers:"));
        summary.addCell(TableUtil.createValueCell(disruptionDriverString()));
        summary.addCell(TableUtil.createKeyCell("Fusion drivers:"));
        summary.addCell(TableUtil.createValueCell(fusionDriverString()));
        summary.addCell(TableUtil.createKeyCell("Integrated viruses:"));
        summary.addCell(TableUtil.createValueCell(virusString()));
        summary.addCell(TableUtil.createKeyCell("Microsatellite indels per Mb:"));
        summary.addCell(TableUtil.createValueCell(msiString()));
        summary.addCell(TableUtil.createKeyCell("Tumor mutational load:"));
        summary.addCell(TableUtil.createValueCell(tmlString()));
        summary.addCell(TableUtil.createKeyCell("CHORD score:"));
        summary.addCell(TableUtil.createValueCell(chordString()));
        summary.addCell(TableUtil.createKeyCell("Tumor mutations per Mb:"));
        summary.addCell(TableUtil.createValueCell(SINGLE_DIGIT.format(report.purple().tumorMutationalBurdenPerMb())));
        summary.addCell(TableUtil.createKeyCell("Number of SVs:"));
        summary.addCell(TableUtil.createValueCell(Integer.toString(report.purple().svTumorMutationalBurden())));
        summary.addCell(TableUtil.createKeyCell("On-label reported treatments:"));
        summary.addCell(TableUtil.createValueCell("todo"));
        summary.addCell(TableUtil.createKeyCell("Off-label reported treatments:"));
        summary.addCell(TableUtil.createValueCell("todo"));

        String circosPath = report.plots().purpleComprehensiveCircosPlot();
        Image circosImage;
        try {
            circosImage = new Image(ImageDataFactory.create(circosPath));
        } catch (MalformedURLException e) {
            throw new IOException("Could not read circos from " + circosPath);
        }
        circosImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
        circosImage.setMaxHeight(300);

        topTable.addCell(summary);
        topTable.addCell(circosImage);

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1 })).setWidth(ReportResources.CONTENT_WIDTH_WIDE).setPadding(0);
        table.addCell(topTable);

        String clonalityPath = report.plots().purpleClonalityPlot();
        Image clonalityImage;
        try {
            clonalityImage = new Image(ImageDataFactory.create(clonalityPath));
        } catch (MalformedURLException e) {
            throw new IOException("Could not read clonality image from " + clonalityPath);
        }
        clonalityImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
        clonalityImage.setMaxHeight(300);

        table.addCell(clonalityImage);
        document.add(table);
    }

    @NotNull
    private String germlineDriverString() {
        return variantDriverString(report.purple().reportableGermlineVariants());
    }

    @NotNull
    private String somaticDriverString() {
        return variantDriverString(report.purple().reportableSomaticVariants());
    }

    @NotNull
    private static String variantDriverString(@NotNull List<ReportableVariant> variants) {
        if (variants.isEmpty()) {
            return "None";
        } else {
            Set<String> highDriverGenes = Sets.newHashSet();
            for (ReportableVariant variant : variants) {
                if (variant.driverLikelihoodInterpretation() == DriverInterpretation.HIGH) {
                    highDriverGenes.add(variant.gene());
                }
            }

            StringJoiner joiner = new StringJoiner(", ");
            for (String gene : highDriverGenes) {
                joiner.add(gene);
            }

            return variants.size() + " (" + joiner.toString() + ")";
        }
    }

    @NotNull
    private String copyNumberDriverString() {
        if (report.purple().reportableGainsLosses().isEmpty()) {
            return "None";
        } else {
            StringJoiner joiner = new StringJoiner(", ");
            for (ReportableGainLoss gainLoss : report.purple().reportableGainsLosses()) {
                joiner.add(gainLoss.gene());
            }
            return report.purple().reportableGainsLosses().size() + " (" + joiner.toString() + ")";
        }
    }

    @NotNull
    private String disruptionDriverString() {
        if (report.linx().homozygousDisruptions().isEmpty()) {
            return "None";
        } else {
            StringJoiner joiner = new StringJoiner(", ");
            for (ReportableHomozygousDisruption disruption : report.linx().homozygousDisruptions()) {
                joiner.add(disruption.gene());
            }
            return report.linx().homozygousDisruptions().size() + " (" + joiner.toString() + ")";
        }
    }

    @NotNull
    private String fusionDriverString() {
        if (report.linx().reportableFusions().isEmpty()) {
            return "None";
        } else {
            StringJoiner joiner = new StringJoiner(", ");
            for (LinxFusion fusion : report.linx().reportableFusions()) {
                joiner.add(fusion.name());
            }
            return report.linx().reportableFusions().size() + " (" + joiner.toString() + ")";
        }
    }

    @NotNull
    private String virusString() {
        if (report.virusInterpreter().reportableViruses().isEmpty()) {
            return "None";
        } else {
            Set<String> viruses = Sets.newHashSet();
            for (AnnotatedVirus virus : report.virusInterpreter().reportableViruses()) {
                if (virus.interpretation() != null) {
                    viruses.add(virus.interpretation().toString());
                } else {
                    viruses.add(virus.name());
                }
            }

            StringJoiner joiner = new StringJoiner(", ");
            for (String virus : viruses) {
                joiner.add(virus);
            }
            return report.virusInterpreter().reportableViruses().size() + " (" + joiner.toString() + ")";
        }
    }

    @NotNull
    private String purpleQCString() {
        StringJoiner joiner = new StringJoiner(", ");
        for (PurpleQCStatus status : report.purple().purpleQC()) {
            joiner.add(status.toString());
        }
        return joiner.toString();
    }

    @NotNull
    private String purityString() {
        DecimalFormat purityFormat = ReportResources.decimalFormat("#'%'");
        return String.format("%s (%s-%s)",
                purityFormat.format(report.purple().purity() * 100),
                purityFormat.format(report.purple().minPurity() * 100),
                purityFormat.format(report.purple().maxPurity() * 100));
    }

    @NotNull
    private String ploidyString() {
        DecimalFormat ploidyFormat = ReportResources.decimalFormat("#.##");
        return String.format("%s (%s-%s)",
                ploidyFormat.format(report.purple().ploidy()),
                ploidyFormat.format(report.purple().minPloidy()),
                ploidyFormat.format(report.purple().maxPloidy()));
    }

    @NotNull
    private String tmlString() {
        return report.purple().tumorMutationalLoad() + " (" + report.purple().tumorMutationalLoadStatus().display() + ")";
    }

    @NotNull
    private String msiString() {
        return SINGLE_DIGIT.format(report.purple().microsatelliteIndelsPerMb()) + " (" + report.purple().microsatelliteStatus().display()
                + ")";
    }

    @NotNull
    private String chordString() {
        return SINGLE_DIGIT.format(report.chord().hrdValue()) + " (" + report.chord().hrStatus().display() + ")";
    }

    @NotNull
    private static String toConfiguredTumorType(@NotNull Set<DoidNode> nodes) {
        StringJoiner joiner = new StringJoiner(", ");

        for (DoidNode node : nodes) {
            joiner.add(node.doidTerm() + " (DOID " + node.doid() + ")");
        }

        return joiner.toString();
    }
}
