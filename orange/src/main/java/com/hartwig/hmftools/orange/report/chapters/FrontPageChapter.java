package com.hartwig.hmftools.orange.report.chapters;

import java.text.DecimalFormat;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Stream;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.datamodel.chord.ChordRecord;
import com.hartwig.hmftools.datamodel.chord.ChordStatus;
import com.hartwig.hmftools.datamodel.cohort.Evaluation;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.cuppa.CuppaPrediction;
import com.hartwig.hmftools.datamodel.linx.HomozygousDisruption;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.orange.OrangeDoidNode;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.orange.PercentileType;
import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.datamodel.purple.PurpleCharacteristics;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleGainLoss;
import com.hartwig.hmftools.datamodel.purple.PurpleMicrosatelliteStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleTumorMutationalStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.virus.AnnotatedVirus;
import com.hartwig.hmftools.datamodel.virus.VirusInterpretation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.orange.algo.cuppa.CuppaInterpretation;
import com.hartwig.hmftools.orange.algo.purple.DriverInterpretation;
import com.hartwig.hmftools.orange.cohort.mapping.CohortConstants;
import com.hartwig.hmftools.orange.report.PlotPathResolver;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Drivers;
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
import org.jetbrains.annotations.Nullable;

public class FrontPageChapter implements ReportChapter {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#.#");
    private static final DecimalFormat TWO_DIGITS = ReportResources.decimalFormat("#.##");
    private static final DecimalFormat PERCENTAGE = ReportResources.decimalFormat("#'%'");

    private static final String NONE = "None";

    @NotNull
    private final OrangeRecord report;
    @NotNull
    private final PlotPathResolver plotPathResolver;
    @NotNull
    private final ReportResources reportResources;

    public FrontPageChapter(@NotNull final OrangeRecord report, @NotNull final PlotPathResolver plotPathResolver,
            @NotNull final ReportResources reportResources) {
        this.report = report;
        this.plotPathResolver = plotPathResolver;
        this.reportResources = reportResources;
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
        Cells cells = new Cells(reportResources);
        Table table = Tables.createContent(contentWidth(),
                new float[] { 3, 2, 1 },
                new Cell[] { cells.createHeader("Configured Primary Tumor"), cells.createHeader("Cuppa Cancer Type"),
                        cells.createHeader("QC") });

        table.addCell(cells.createContent(configuredPrimaryTumor(report.configuredPrimaryTumor())));
        table.addCell(cells.createContent(cuppaCancerType(report.cuppa())));
        table.addCell(cells.createContent(purpleQCString()));
        document.add(new Tables(reportResources).createWrapping(table));
    }

    @NotNull
    private static String cuppaCancerType(@Nullable CuppaData cuppa) {
        if (cuppa == null) {
            return ReportResources.NOT_AVAILABLE;
        }

        CuppaPrediction best = CuppaInterpretation.best(cuppa);
        return best.cancerType() + " (" + PERCENTAGE.format(best.likelihood() * 100) + ")";
    }

    @NotNull
    private static String configuredPrimaryTumor(@NotNull Set<OrangeDoidNode> nodes) {
        Set<String> configured = Sets.newHashSet();
        for (OrangeDoidNode node : nodes) {
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
        Cells cells = new Cells(reportResources);
        Stream.of(Maps.immutableEntry("Purity:", purityString()),
                Maps.immutableEntry("Ploidy:", ploidyString()),
                Maps.immutableEntry("Somatic variant drivers:", somaticVariantDriverString()),
                Maps.immutableEntry("Germline variant drivers:", germlineVariantDriverString()),
                Maps.immutableEntry("Somatic copy number drivers:", somaticCopyNumberDriverString()),
                Maps.immutableEntry("Germline copy number drivers:", germlineCopyNumberDriverString()),
                Maps.immutableEntry("Somatic disruption drivers:", somaticDisruptionDriverString()),
                Maps.immutableEntry("Germline disruption drivers:", germlineDisruptionDriverString()),
                Maps.immutableEntry("Fusion drivers:", fusionDriverString()),
                Maps.immutableEntry("Viral presence:", virusString()),
                Maps.immutableEntry("Whole genome duplicated:", report.purple().characteristics().wholeGenomeDuplication() ? "Yes" : "No"),
                Maps.immutableEntry("Microsatellite indels per Mb:", msiString()),
                Maps.immutableEntry("Tumor mutations per Mb:",
                        SINGLE_DIGIT.format(report.purple().characteristics().tumorMutationalBurdenPerMb())),
                Maps.immutableEntry("Tumor mutational load:", tmlString()),
                Maps.immutableEntry("HR deficiency score:", hrDeficiencyString()),
                Maps.immutableEntry("DPYD status:", dpydStatus()),
                Maps.immutableEntry("Number of SVs:", svTmbString()),
                Maps.immutableEntry("Max complex cluster size:", maxComplexSizeString()),
                Maps.immutableEntry("Telomeric SGLs:", telomericSGLString()),
                Maps.immutableEntry("Number of LINE insertions:", lineCountString())).forEach(entry -> {
            summary.addCell(cells.createKey(entry.getKey()));
            summary.addCell(cells.createValue(entry.getValue()));
        });

        Image circosImage = Images.build(plotPathResolver.resolve(report.plots().purpleFinalCircosPlot()));
        circosImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
        circosImage.setMaxHeight(290);

        topTable.addCell(summary);
        topTable.addCell(circosImage);

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1 })).setWidth(contentWidth()).setPadding(0);
        table.addCell(topTable);

        Image clonalityImage = Images.build(plotPathResolver.resolve(report.plots().purpleClonalityPlot()));
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
    private String somaticVariantDriverString() {
        return variantDriverString(report.purple().reportableSomaticVariants(), report.purple().somaticDrivers());
    }

    @NotNull
    private String germlineVariantDriverString() {
        List<PurpleVariant> reportableGermlineVariants = report.purple().reportableGermlineVariants();
        List<PurpleDriver> germlineDrivers = report.purple().germlineDrivers();
        if (reportableGermlineVariants != null && germlineDrivers != null) {
            return variantDriverString(reportableGermlineVariants, germlineDrivers);
        } else {
            return ReportResources.NOT_AVAILABLE;
        }
    }

    @NotNull
    private static String variantDriverString(@NotNull List<PurpleVariant> variants, @NotNull List<PurpleDriver> drivers) {
        if (variants.isEmpty()) {
            return NONE;
        }

        Set<String> highDriverGenes = Sets.newTreeSet(Comparator.naturalOrder());
        for (PurpleVariant variant : variants) {
            PurpleDriver driver = Drivers.canonicalMutationEntryForGene(drivers, variant.gene());
            if (driver != null && DriverInterpretation.interpret(driver.driverLikelihood()) == DriverInterpretation.HIGH) {
                highDriverGenes.add(variant.gene());
            }
        }

        return !highDriverGenes.isEmpty() ? variants.size() + " (" + concat(highDriverGenes) + ")" : String.valueOf(variants.size());
    }

    @NotNull
    private String somaticCopyNumberDriverString() {
        return copyNumberDriverString(report.purple().reportableSomaticGainsLosses());
    }

    @NotNull
    private String germlineCopyNumberDriverString() {
        List<PurpleGainLoss> germlineGainsLosses = report.purple().reportableGermlineFullLosses();
        if (germlineGainsLosses == null) {
            return ReportResources.NOT_AVAILABLE;
        }
        return copyNumberDriverString(germlineGainsLosses);
    }

    @NotNull
    private static String copyNumberDriverString(@NotNull List<PurpleGainLoss> gainsLosses) {
        if (gainsLosses.isEmpty()) {
            return NONE;
        }

        Set<String> genes = Sets.newTreeSet(Comparator.naturalOrder());
        for (PurpleGainLoss gainLoss : gainsLosses) {
            genes.add(gainLoss.gene());
        }
        return gainsLosses.size() + " (" + concat(genes) + ")";
    }

    @NotNull
    private String somaticDisruptionDriverString() {
        return disruptionDriverString(report.linx().somaticHomozygousDisruptions());
    }

    @NotNull
    private String germlineDisruptionDriverString() {
        List<HomozygousDisruption> germlineHomozygousDisruptions = report.linx().germlineHomozygousDisruptions();
        if (germlineHomozygousDisruptions == null) {
            return ReportResources.NOT_AVAILABLE;
        }
        return disruptionDriverString(germlineHomozygousDisruptions);
    }

    @NotNull
    private static String disruptionDriverString(@NotNull List<HomozygousDisruption> homozygousDisruptions) {
        if (homozygousDisruptions.isEmpty()) {
            return NONE;
        }

        Set<String> genes = Sets.newTreeSet(Comparator.naturalOrder());
        for (HomozygousDisruption homozygousDisruption : homozygousDisruptions) {
            genes.add(homozygousDisruption.gene());
        }
        return homozygousDisruptions.size() + " (" + concat(genes) + ")";
    }

    @NotNull
    private String fusionDriverString() {
        if (report.linx().reportableSomaticFusions().isEmpty()) {
            return NONE;
        }

        Set<String> fusions = Sets.newTreeSet(Comparator.naturalOrder());
        for (LinxFusion fusion : report.linx().reportableSomaticFusions()) {
            fusions.add(fusion.name());
        }
        return report.linx().reportableSomaticFusions().size() + " (" + concat(fusions) + ")";
    }

    @NotNull
    private String virusString() {
        VirusInterpreterData virusInterpreter = report.virusInterpreter();
        if (virusInterpreter == null) {
            return ReportResources.NOT_AVAILABLE;
        }

        if (virusInterpreter.reportableViruses().isEmpty()) {
            return NONE;
        }

        Set<String> viruses = Sets.newTreeSet(Comparator.naturalOrder());
        for (AnnotatedVirus virus : virusInterpreter.reportableViruses()) {
            VirusInterpretation interpretation = virus.interpretation();
            if (interpretation != null) {
                viruses.add(interpretation.name());
            } else {
                viruses.add(virus.name());
            }
        }

        return virusInterpreter.reportableViruses().size() + " (" + concat(viruses) + ")";
    }

    @NotNull
    private String msiString() {
        PurpleCharacteristics characteristics = report.purple().characteristics();
        return SINGLE_DIGIT.format(characteristics.microsatelliteIndelsPerMb()) + " (" + display(characteristics.microsatelliteStatus()) + ")";
    }

    @NotNull
    private String tmlString() {
        PurpleCharacteristics characteristics = report.purple().characteristics();
        return characteristics.tumorMutationalLoad() + " (" + display(characteristics.tumorMutationalLoadStatus()) + ")";
    }

    private static String display(PurpleMicrosatelliteStatus microsatelliteStatus) {
        switch (microsatelliteStatus) {
            case MSI:
                return "Unstable";
            case MSS:
                return "Stable";
            case UNKNOWN:
                return "Unknown";
        }
        throw new IllegalStateException();
    }

    private static String display(@NotNull PurpleTumorMutationalStatus tumorMutationalStatus) {
        switch (tumorMutationalStatus)
        {
            case HIGH:
                return "High";
            case LOW:
                return "Low";
            case UNKNOWN:
                return "Unknown";
        }
        throw new IllegalStateException();
    }

    @NotNull
    private String hrDeficiencyString() {
        ChordRecord chord = report.chord();
        if (chord == null) {
            return ReportResources.NOT_AVAILABLE;
        }

        if (chord.hrStatus() == ChordStatus.CANNOT_BE_DETERMINED) {
            return displayChordStatus(ChordStatus.CANNOT_BE_DETERMINED);
        }

        String addon = Strings.EMPTY;
        if (chord.hrStatus() == ChordStatus.HR_DEFICIENT) {
            if (chord.hrdType().contains("BRCA1")) {
                addon = " - BRCA1 (" + TWO_DIGITS.format(chord.brca1Value()) + ")";
            } else if (chord.hrdType().contains("BRCA2")) {
                addon = " - BRCA2 (" + TWO_DIGITS.format(chord.brca2Value()) + ")";
            } else {
                addon = chord.hrdType();
            }
        }
        return SINGLE_DIGIT.format(chord.hrdValue()) + " (" + displayChordStatus(chord.hrStatus()) + addon + ")";
    }

    private static String displayChordStatus(ChordStatus chordStatus) {
        switch (chordStatus) {
            case CANNOT_BE_DETERMINED: return "Cannot be determined";
            case HR_PROFICIENT: return "Proficient";
            case HR_DEFICIENT: return "Deficient";
            default: return "Unknown";
        }
    }

    @NotNull
    private String dpydStatus() {
        Set<PeachGenotype> genotypes = report.peach();
        if (genotypes == null) {
            return ReportResources.NOT_AVAILABLE;
        }

        Set<String> haplotypes = Sets.newHashSet();
        for (PeachGenotype genotype : genotypes) {
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
                String cancerTypePercentile = percentile != null ? PERCENTAGE.format(percentile * 100) : "NA";
                addon = addon + " | " + evaluation.cancerType() + " " + cancerTypePercentile;
            }
            addon = addon + ")";
        }

        return svTmb + addon;
    }

    @NotNull
    private String maxComplexSizeString() {
        CuppaData cuppa = report.cuppa();
        return cuppa != null ? Integer.toString(cuppa.maxComplexSize()) : ReportResources.NOT_AVAILABLE;
    }

    @NotNull
    private String telomericSGLString() {
        CuppaData cuppa = report.cuppa();
        return cuppa != null ? Integer.toString(cuppa.telomericSGLs()) : ReportResources.NOT_AVAILABLE;
    }

    @NotNull
    private String lineCountString() {
        CuppaData cuppa = report.cuppa();
        return cuppa != null ? Integer.toString(cuppa.lineCount()) : ReportResources.NOT_AVAILABLE;
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
