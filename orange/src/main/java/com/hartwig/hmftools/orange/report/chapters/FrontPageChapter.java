package com.hartwig.hmftools.orange.report.chapters;

import java.text.DecimalFormat;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.common.linx.HomozygousDisruption;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.interpretation.GainLoss;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.DriverInterpretation;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantFactory;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.algo.cuppa.CuppaData;
import com.hartwig.hmftools.orange.algo.cuppa.CuppaInterpretation;
import com.hartwig.hmftools.orange.algo.cuppa.CuppaPrediction;
import com.hartwig.hmftools.orange.algo.protect.ProtectInterpretedData;
import com.hartwig.hmftools.orange.algo.purple.PurpleCharacteristics;
import com.hartwig.hmftools.orange.cohort.datamodel.Evaluation;
import com.hartwig.hmftools.orange.cohort.mapping.CohortConstants;
import com.hartwig.hmftools.orange.cohort.percentile.PercentileType;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Images;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.HorizontalAlignment;
import com.itextpdf.layout.property.UnitValue;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class FrontPageChapter implements ReportChapter {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#.#");
    private static final DecimalFormat TWO_DIGITS = ReportResources.decimalFormat("#.##");
    private static final DecimalFormat PERCENTAGE = ReportResources.decimalFormat("#'%'");

    private static final String NONE = "None";

    @NotNull
    private final OrangeReport report;
    private final boolean reportGermline;

    public FrontPageChapter(@NotNull final OrangeReport report, final boolean reportGermline) {
        this.report = report;
        this.reportGermline = reportGermline;
    }

    @NotNull
    @Override
    public String name() {
        return "Front Page";
    }

    @NotNull
    @Override
    public PageSize pageSize() {
        return PageSize.A4;
    }

    @Override
    public void render(@NotNull Document document) {
        addSummaryTable(document);
        addDetailsAndPlots(document);
    }

    private void addSummaryTable(@NotNull Document document) {
        Table table = Tables.createContent(contentWidth(),
                new float[] { 3, 2, 1 },
                new Cell[] { Cells.createHeader("Configured Primary Tumor"), Cells.createHeader("Cuppa Cancer Type"),
                        Cells.createHeader("QC") });

        table.addCell(Cells.createContent(configuredPrimaryTumor(report.configuredPrimaryTumor())));
        table.addCell(Cells.createContent(cuppaCancerType(report.cuppa())));
        table.addCell(Cells.createContent(purpleQCString()));
        document.add(Tables.createWrapping(table));
    }

    @NotNull
    private static String cuppaCancerType(@NotNull CuppaData cuppa) {
        CuppaPrediction best = CuppaInterpretation.best(cuppa);
        return best.cancerType() + " (" + PERCENTAGE.format(best.likelihood() * 100) + ")";
    }

    @NotNull
    private static String configuredPrimaryTumor(@NotNull Set<DoidNode> nodes) {
        Set<String> configured = Sets.newHashSet();
        for (DoidNode node : nodes) {
            configured.add(node.doidTerm() + " (DOID " + node.doid() + ")");
        }

        return concat(configured);
    }

    @NotNull
    private String purpleQCString() {
        Set<String> purpleStatuses = Sets.newHashSet();
        for (PurpleQCStatus status : report.purple().fit().qc().status()) {
            purpleStatuses.add(status.toString());
        }
        return concat(purpleStatuses);
    }

    private void addDetailsAndPlots(@NotNull Document document) {
        Table topTable = new Table(UnitValue.createPercentArray(new float[] { 1, 1 })).setWidth(contentWidth() - 5);

        Table summary = new Table(UnitValue.createPercentArray(new float[] { 1, 1 }));
        summary.addCell(Cells.createKey("Purity:"));
        summary.addCell(Cells.createValue(purityString()));
        summary.addCell(Cells.createKey("Ploidy:"));
        summary.addCell(Cells.createValue(ploidyString()));
        summary.addCell(Cells.createKey("Somatic variant drivers:"));
        summary.addCell(Cells.createValue(somaticDriverString()));
        summary.addCell(Cells.createKey("Germline variant drivers:"));
        summary.addCell(Cells.createValue(germlineDriverString()));
        summary.addCell(Cells.createKey("Copy number drivers:"));
        summary.addCell(Cells.createValue(copyNumberDriverString()));
        summary.addCell(Cells.createKey("Disruption drivers:"));
        summary.addCell(Cells.createValue(disruptionDriverString()));
        summary.addCell(Cells.createKey("Fusion drivers:"));
        summary.addCell(Cells.createValue(fusionDriverString()));
        summary.addCell(Cells.createKey("Viral presence:"));
        summary.addCell(Cells.createValue(virusString()));
        summary.addCell(Cells.createKey("Whole genome duplicated:"));
        summary.addCell(Cells.createValue(report.purple().characteristics().wholeGenomeDuplication() ? "Yes" : "No"));
        summary.addCell(Cells.createKey("Microsatellite indels per Mb:"));
        summary.addCell(Cells.createValue(msiString()));
        summary.addCell(Cells.createKey("Tumor mutations per Mb:"));
        summary.addCell(Cells.createValue(SINGLE_DIGIT.format(report.purple().characteristics().tumorMutationalBurdenPerMb())));
        summary.addCell(Cells.createKey("Tumor mutational load:"));
        summary.addCell(Cells.createValue(tmlString()));
        summary.addCell(Cells.createKey("HR deficiency score:"));
        summary.addCell(Cells.createValue(hrDeficiencyString()));
        summary.addCell(Cells.createKey("DPYD status:"));
        summary.addCell(Cells.createValue(dpydStatus()));
        summary.addCell(Cells.createKey("Number of SVs:"));
        summary.addCell(Cells.createValue(svTmbString()));
        summary.addCell(Cells.createKey("Max complex cluster size:"));
        summary.addCell(Cells.createValue(Integer.toString(report.cuppa().maxComplexSize())));
        summary.addCell(Cells.createKey("Telomeric SGLs:"));
        summary.addCell(Cells.createValue(Integer.toString(report.cuppa().telomericSGLs())));
        summary.addCell(Cells.createKey("Number of LINE insertions:"));
        summary.addCell(Cells.createValue(Integer.toString(report.cuppa().LINECount())));
        summary.addCell(Cells.createKey("On-label treatments:"));
        summary.addCell(Cells.createValue(onLabelTreatmentString()));
        summary.addCell(Cells.createKey("Off-label treatments:"));
        summary.addCell(Cells.createValue(offLabelTreatmentString()));

        Image circosImage = Images.build(report.plots().purpleFinalCircosPlot());
        circosImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
        circosImage.setMaxHeight(290);

        topTable.addCell(summary);
        topTable.addCell(circosImage);

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1 })).setWidth(contentWidth()).setPadding(0);
        table.addCell(topTable);

        Image clonalityImage = Images.build(report.plots().purpleClonalityPlot());
        clonalityImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
        clonalityImage.setMaxHeight(270);

        table.addCell(clonalityImage);
        document.add(table);
    }

    @NotNull
    private String purityString() {
        return String.format("%s (%s-%s)",
                PERCENTAGE.format(report.purple().fit().purity() * 100),
                PERCENTAGE.format(report.purple().fit().minPurity() * 100),
                PERCENTAGE.format(report.purple().fit().maxPurity() * 100));
    }

    @NotNull
    private String ploidyString() {
        return String.format("%s (%s-%s)",
                TWO_DIGITS.format(report.purple().fit().ploidy()),
                TWO_DIGITS.format(report.purple().fit().minPloidy()),
                TWO_DIGITS.format(report.purple().fit().maxPloidy()));
    }

    @NotNull
    private String somaticDriverString() {
        List<ReportableVariant> reportableVariants;
        if (reportGermline) {
            reportableVariants = report.purple().reportableSomaticVariants();
        } else {
            reportableVariants = ReportableVariantFactory.mergeVariantLists(report.purple().reportableSomaticVariants(),
                    report.purple().reportableGermlineVariants());
        }
        return variantDriverString(reportableVariants);
    }

    @NotNull
    private String germlineDriverString() {
        if (reportGermline) {
            return variantDriverString(report.purple().reportableGermlineVariants());
        } else {
            return ReportResources.NOT_AVAILABLE;
        }
    }

    @NotNull
    private static String variantDriverString(@NotNull List<ReportableVariant> variants) {
        if (variants.isEmpty()) {
            return NONE;
        }

        Set<String> highDriverGenes = Sets.newTreeSet(Comparator.naturalOrder());
        for (ReportableVariant variant : variants) {
            if (variant.driverLikelihoodInterpretation() == DriverInterpretation.HIGH) {
                highDriverGenes.add(variant.gene());
            }
        }

        return !highDriverGenes.isEmpty() ? variants.size() + " (" + concat(highDriverGenes) + ")" : String.valueOf(variants.size());
    }

    @NotNull
    private String copyNumberDriverString() {
        if (report.purple().reportableSomaticGainsLosses().isEmpty()) {
            return NONE;
        }

        Set<String> genes = Sets.newTreeSet(Comparator.naturalOrder());
        for (GainLoss gainLoss : report.purple().reportableSomaticGainsLosses()) {
            genes.add(gainLoss.gene());
        }
        return report.purple().reportableSomaticGainsLosses().size() + " (" + concat(genes) + ")";
    }

    @NotNull
    private String disruptionDriverString() {
        if (report.linx().homozygousDisruptions().isEmpty()) {
            return NONE;
        }

        Set<String> genes = Sets.newTreeSet(Comparator.naturalOrder());
        for (HomozygousDisruption disruption : report.linx().homozygousDisruptions()) {
            genes.add(disruption.gene());
        }
        return report.linx().homozygousDisruptions().size() + " (" + concat(genes) + ")";
    }

    @NotNull
    private String fusionDriverString() {
        if (report.linx().reportableFusions().isEmpty()) {
            return NONE;
        }

        Set<String> fusions = Sets.newTreeSet(Comparator.naturalOrder());
        for (LinxFusion fusion : report.linx().reportableFusions()) {
            fusions.add(fusion.name());
        }
        return report.linx().reportableFusions().size() + " (" + concat(fusions) + ")";
    }

    @NotNull
    private String virusString() {
        if (report.virusInterpreter().reportableViruses().isEmpty()) {
            return NONE;
        }

        Set<String> viruses = Sets.newTreeSet(Comparator.naturalOrder());
        for (AnnotatedVirus virus : report.virusInterpreter().reportableViruses()) {
            if (virus.interpretation() != null) {
                viruses.add(virus.interpretation());
            } else {
                viruses.add(virus.name());
            }
        }

        return report.virusInterpreter().reportableViruses().size() + " (" + concat(viruses) + ")";
    }

    @NotNull
    private String msiString() {
        PurpleCharacteristics characteristics = report.purple().characteristics();
        return SINGLE_DIGIT.format(characteristics.microsatelliteIndelsPerMb()) + " (" + characteristics.microsatelliteStatus().display()
                + ")";
    }

    @NotNull
    private String tmlString() {
        PurpleCharacteristics characteristics = report.purple().characteristics();
        return characteristics.tumorMutationalLoad() + " (" + characteristics.tumorMutationalLoadStatus().display() + ")";
    }

    @NotNull
    private String hrDeficiencyString() {
        ChordData chord = report.chord();
        if (chord.hrStatus() == ChordStatus.CANNOT_BE_DETERMINED) {
            return ChordStatus.CANNOT_BE_DETERMINED.display();
        }

        String addon = Strings.EMPTY;
        if (chord.hrStatus() == ChordStatus.HR_DEFICIENT) {
            if (chord.hrdType().contains("BRCA1")) {
                addon = " - BRCA1 (" + TWO_DIGITS.format(chord.BRCA1Value()) + ")";
            } else if (chord.hrdType().contains("BRCA2")) {
                addon = " - BRCA2 (" + TWO_DIGITS.format(chord.BRCA2Value()) + ")";
            } else {
                addon = chord.hrdType();
            }
        }
        return SINGLE_DIGIT.format(chord.hrdValue()) + " (" + chord.hrStatus().display() + addon + ")";
    }

    @NotNull
    private String dpydStatus() {
        Set<String> haplotypes = Sets.newHashSet();
        for (PeachGenotype genotype : report.peach()) {
            if (genotype.gene().equals("DPYD")) {
                haplotypes.add(genotype.haplotype() + " (" + genotype.function() + ")");
            }
        }
        return !haplotypes.isEmpty() ? concat(haplotypes) : NONE;
    }

    @NotNull
    private String svTmbString() {
        String svTmb = String.valueOf(report.purple().characteristics().svTumorMutationalBurden());

        Evaluation evaluation = report.cohortEvaluations().get(PercentileType.SV_TMB);
        String addon = Strings.EMPTY;
        if (evaluation != null) {
            String panCancerPercentile = PERCENTAGE.format(evaluation.panCancerPercentile() * 100);
            addon = " (Pan " + panCancerPercentile;
            String cancerType = evaluation.cancerType();
            if (cancerType != null && !cancerType.equals(CohortConstants.COHORT_OTHER)
                    && !cancerType.equals(CohortConstants.COHORT_UNKNOWN)) {
                Double percentile = evaluation.cancerTypePercentile();
                String cancerTypePercentile = percentile != null ? PERCENTAGE.format(evaluation.cancerTypePercentile() * 100) : "NA";
                addon = addon + " | " + evaluation.cancerType() + " " + cancerTypePercentile;
            }
            addon = addon + ")";
        }

        return svTmb + addon;
    }

    @NotNull
    private String onLabelTreatmentString() {
        return treatmentString(report.protect(), true);
    }

    @NotNull
    private String offLabelTreatmentString() {
        return treatmentString(report.protect(), false);
    }

    @NotNull
    private static String treatmentString(@NotNull ProtectInterpretedData protect, boolean requireOnLabel) {
        Set<String> levels = Sets.newTreeSet(Comparator.naturalOrder());
        Set<String> treatments = Sets.newHashSet();
        for (ProtectEvidence evidence : protect.reportableEvidences()) {
            if (evidence.onLabel() == requireOnLabel) {
                treatments.add(evidence.treatment());
                levels.add(evidence.level().toString());
            }
        }

        for (ProtectEvidence evidence : protect.unreportedEvidences()) {
            if (evidence.onLabel() == requireOnLabel) {
                treatments.add(evidence.treatment());
                levels.add(evidence.level().toString());
            }
        }

        if (treatments.isEmpty()) {
            return NONE;
        } else {
            return treatments.size() + " (" + concat(levels) + ")";
        }
    }

    @NotNull
    private static String concat(@NotNull Iterable<String> strings) {
        StringJoiner joiner = new StringJoiner(", ");
        for (String string : strings) {
            joiner.add(string);
        }

        return joiner.toString();
    }
}
