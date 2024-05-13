package com.hartwig.hmftools.orange.report.chapters;

import static com.hartwig.hmftools.common.peach.PeachUtil.UNKNOWN_ALLELE_STRING;
import static com.hartwig.hmftools.common.peach.PeachUtil.convertToZygosityString;
import static com.hartwig.hmftools.orange.report.ReportResources.formatPercentage;
import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.ReportResources.formatTwoDigitDecimal;

import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.datamodel.chord.ChordRecord;
import com.hartwig.hmftools.datamodel.chord.ChordStatus;
import com.hartwig.hmftools.datamodel.cohort.Evaluation;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.cuppa.CuppaPrediction;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.linx.LinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.orange.OrangeDoidNode;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.orange.PercentileType;
import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.datamodel.purple.PurpleCharacteristics;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleGainLoss;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.virus.VirusInterpretation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;
import com.hartwig.hmftools.orange.algo.purple.DriverInterpretation;
import com.hartwig.hmftools.orange.cohort.mapping.CohortConstants;
import com.hartwig.hmftools.orange.report.PlotPathResolver;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Drivers;
import com.hartwig.hmftools.orange.report.interpretation.PurpleQCInterpretation;
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

public class FrontPageChapter implements ReportChapter
{
    private static final String NONE = "None";

    @NotNull
    private final OrangeRecord report;
    @NotNull
    private final PlotPathResolver plotPathResolver;
    @NotNull
    private final ReportResources reportResources;

    public FrontPageChapter(@NotNull final OrangeRecord report, @NotNull final PlotPathResolver plotPathResolver,
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
        return "Front Page";
    }

    @NotNull
    @Override
    public PageSize pageSize()
    {
        return PageSize.A4;
    }

    @Override
    public void render(@NotNull Document document)
    {
        addSummaryTable(document);
        addDetailsAndPlots(document);
    }

    private void addSummaryTable(@NotNull Document document)
    {
        Cells cells = new Cells(reportResources);

        float[] headerComponents = new float[] { 3, 2, 1 };

        Cell[] headerCells = new Cell[] {
                cells.createHeader("Configured Primary Tumor"),
                cells.createHeader(!report.tumorOnlyMode() ? "Cuppa Cancer Type" : "Tumor-Only"),
                cells.createHeader("QC") };

        Table table = Tables.createContent(contentWidth(), headerComponents, headerCells);

        table.addCell(cells.createContent(configuredPrimaryTumor(report.configuredPrimaryTumor())));
        table.addCell(cells.createContent(cuppaCancerTypeString()));
        table.addCell(cells.createContent(purpleQCString()));

        addQCWarningInCaseOfFail(table, cells);

        document.add(new Tables(reportResources).createWrapping(table));
    }

    private void addQCWarningInCaseOfFail(@NotNull Table table, @NotNull Cells cells)
    {
        boolean isFailNoTumor = PurpleQCInterpretation.isFailNoTumor(report.purple().fit().qc());
        boolean isContaminated = PurpleQCInterpretation.isContaminated(report.purple().fit().qc());

        if(isFailNoTumor || isContaminated)
        {
            String reason;
            if(isFailNoTumor && isContaminated)
            {
                reason = "no tumor and contamination";
            }
            else if(isFailNoTumor)
            {
                reason = "no tumor";
            }
            else
            {
                reason = "contamination";
            }

            String warning = "The QC status of this sample is fail (" + reason + ")" +
                    ": all presented data in this report should be interpreted with caution";
            table.addCell(cells.createSpanningWarning(table, warning));
        }
    }

    @NotNull
    private static String configuredPrimaryTumor(@NotNull Set<OrangeDoidNode> nodes)
    {
        Set<String> configured = Sets.newHashSet();
        for(OrangeDoidNode node : nodes)
        {
            configured.add(node.doidTerm() + " (DOID " + node.doid() + ")");
        }

        return concat(configured);
    }

    private String cuppaCancerTypeString()
    {
        if(report.tumorOnlyMode())
        {
            return "";
        }

        if(PurpleQCInterpretation.isFail(report.purple().fit().qc()))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        return cuppaCancerType(report.cuppa());
    }

    @NotNull
    private static String cuppaCancerType(@Nullable CuppaData cuppaData)
    {
        if(cuppaData == null)
        {
            return ReportResources.NOT_AVAILABLE;
        }

        CuppaPrediction best = cuppaData.bestPrediction();
        return best.cancerType() + " (" + formatPercentage(best.likelihood()) + ")";
    }

    @NotNull
    private String purpleQCString()
    {
        Set<String> purpleStatuses = Sets.newHashSet();
        for(PurpleQCStatus status : report.purple().fit().qc().status())
        {
            purpleStatuses.add(status.toString());
        }
        return concat(purpleStatuses);
    }

    private void addDetailsAndPlots(@NotNull Document document)
    {
        Table topTable = new Table(UnitValue.createPercentArray(new float[] { 1, 1 })).setWidth(contentWidth() - 5);

        Table summary = new Table(UnitValue.createPercentArray(new float[] { 1, 1 }));
        Cells cells = new Cells(reportResources);

        boolean includeGermline = !report.tumorOnlyMode();

        addCellEntry(summary, cells, "Purity:", purityString());
        addCellEntry(summary, cells, "Ploidy:", ploidyString());
        addCellEntry(summary, cells, "Somatic variant drivers:", somaticVariantDriverString());

        if(includeGermline)
        {
            addCellEntry(summary, cells, "Germline variant drivers:", germlineVariantDriverString());
        }

        addCellEntry(summary, cells, "Somatic copy number drivers:", somaticCopyNumberDriverString());

        if(includeGermline)
        {
            addCellEntry(summary, cells, "Germline copy number drivers:", germlineCopyNumberDriverString());
        }

        addCellEntry(summary, cells, "Somatic disruption drivers:", somaticDisruptionDriverString());

        if(includeGermline)
        {
            addCellEntry(summary, cells, "Germline disruption drivers:", germlineDisruptionDriverString());
        }

        addCellEntry(summary, cells, "Fusion drivers:", fusionDriverString());

        if(includeGermline)
        {
            addCellEntry(summary, cells, "Viral presence:", virusString());
        }

        addCellEntry(summary, cells, "Whole genome duplicated:", wgdString());
        addCellEntry(summary, cells, "Microsatellite indels per Mb:", msiString());
        addCellEntry(summary, cells, "Tumor mutations per Mb:", tmbString());
        addCellEntry(summary, cells, "Tumor mutational load:", tmlString());

        if(includeGermline) // will change once we solve HRD for targeted mode
        {
            addCellEntry(summary, cells, "HR deficiency score:", hrDeficiencyString());

            addCellEntry(summary, cells, "DPYD status:", geneStatus("DPYD"));
            addCellEntry(summary, cells, "UGT1A1 status:", geneStatus("UGT1A1"));

            addCellEntry(summary, cells, "Number of SVs:", svTmbString());
            addCellEntry(summary, cells, "Max complex cluster size:", maxComplexSizeString());
            addCellEntry(summary, cells, "Telomeric SGLs:", telomericSGLString());
            addCellEntry(summary, cells, "Number of LINE insertions:", lineCountString());
        }

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

    private static void addCellEntry(final Table summary, final Cells cells, final String name, final String value)
    {
        summary.addCell(cells.createKey(name));
        summary.addCell(cells.createValue(value));
    }

    @NotNull
    private String purityString()
    {
        return String.format("%s (%s-%s)",
                formatPercentage(report.purple().fit().purity()),
                formatPercentage(report.purple().fit().minPurity()),
                formatPercentage(report.purple().fit().maxPurity()));
    }

    @NotNull
    private String ploidyString()
    {
        return String.format("%s (%s-%s)",
                formatTwoDigitDecimal(report.purple().fit().ploidy()),
                formatTwoDigitDecimal(report.purple().fit().minPloidy()),
                formatTwoDigitDecimal(report.purple().fit().maxPloidy()));
    }

    @NotNull
    private String somaticVariantDriverString()
    {
        if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        return variantDriverString(report.purple().reportableSomaticVariants(), report.purple().somaticDrivers());
    }

    @NotNull
    private String germlineVariantDriverString()
    {
        if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        List<PurpleVariant> reportableGermlineVariants = report.purple().reportableGermlineVariants();
        List<PurpleDriver> germlineDrivers = report.purple().germlineDrivers();
        if(reportableGermlineVariants != null && germlineDrivers != null)
        {
            return variantDriverString(reportableGermlineVariants, germlineDrivers);
        }
        else
        {
            return ReportResources.NOT_AVAILABLE;
        }
    }

    @NotNull
    private static String variantDriverString(@NotNull List<PurpleVariant> variants, @NotNull List<PurpleDriver> drivers)
    {
        if(variants.isEmpty())
        {
            return NONE;
        }

        Set<String> highDriverGenes = Sets.newTreeSet(Comparator.naturalOrder());
        for(PurpleVariant variant : variants)
        {
            PurpleDriver driver = Drivers.canonicalMutationEntryForGene(drivers, variant.gene());
            if(driver != null && DriverInterpretation.interpret(driver.driverLikelihood()) == DriverInterpretation.HIGH)
            {
                highDriverGenes.add(variant.gene());
            }
        }

        return !highDriverGenes.isEmpty() ? variants.size() + " (" + concat(highDriverGenes) + ")" : String.valueOf(variants.size());
    }

    @NotNull
    private String somaticCopyNumberDriverString()
    {
        if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        return copyNumberDriverString(report.purple().reportableSomaticGainsLosses());
    }

    @NotNull
    private String germlineCopyNumberDriverString()
    {
        if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        List<PurpleGainLoss> germlineGainsLosses = report.purple().reportableGermlineFullLosses();
        if(germlineGainsLosses == null)
        {
            return ReportResources.NOT_AVAILABLE;
        }
        return copyNumberDriverString(germlineGainsLosses);
    }

    @NotNull
    private static String copyNumberDriverString(@NotNull List<PurpleGainLoss> gainsLosses)
    {
        if(gainsLosses.isEmpty())
        {
            return NONE;
        }

        Set<String> genes = Sets.newTreeSet(Comparator.naturalOrder());
        for(PurpleGainLoss gainLoss : gainsLosses)
        {
            genes.add(gainLoss.gene());
        }
        return gainsLosses.size() + " (" + concat(genes) + ")";
    }

    @NotNull
    private String somaticDisruptionDriverString()
    {
        if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        return disruptionDriverString(report.linx().somaticHomozygousDisruptions());
    }

    @NotNull
    private String germlineDisruptionDriverString()
    {
        if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        List<LinxHomozygousDisruption> germlineHomozygousDisruptions = report.linx().germlineHomozygousDisruptions();
        if(germlineHomozygousDisruptions == null)
        {
            return ReportResources.NOT_AVAILABLE;
        }
        return disruptionDriverString(germlineHomozygousDisruptions);
    }

    @NotNull
    private static String disruptionDriverString(@NotNull List<LinxHomozygousDisruption> homozygousDisruptions)
    {
        if(homozygousDisruptions.isEmpty())
        {
            return NONE;
        }

        Set<String> genes = Sets.newTreeSet(Comparator.naturalOrder());
        for(LinxHomozygousDisruption homozygousDisruption : homozygousDisruptions)
        {
            genes.add(homozygousDisruption.gene());
        }
        return homozygousDisruptions.size() + " (" + concat(genes) + ")";
    }

    @NotNull
    private String fusionDriverString()
    {
        if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        if(report.linx().reportableSomaticFusions().isEmpty())
        {
            return NONE;
        }

        Set<String> fusions = Sets.newTreeSet(Comparator.naturalOrder());
        for(LinxFusion fusion : report.linx().reportableSomaticFusions())
        {
            fusions.add(fusion.display());
        }
        return report.linx().reportableSomaticFusions().size() + " (" + concat(fusions) + ")";
    }

    @NotNull
    private String virusString()
    {
        if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        VirusInterpreterData virusInterpreter = report.virusInterpreter();
        if(virusInterpreter == null)
        {
            return ReportResources.NOT_AVAILABLE;
        }

        if(virusInterpreter.reportableViruses().isEmpty())
        {
            return NONE;
        }

        Set<String> viruses = Sets.newTreeSet(Comparator.naturalOrder());
        for(VirusInterpreterEntry virus : virusInterpreter.reportableViruses())
        {
            VirusInterpretation interpretation = virus.interpretation();
            if(interpretation != null)
            {
                viruses.add(interpretation.name());
            }
            else
            {
                viruses.add(virus.name());
            }
        }

        return virusInterpreter.reportableViruses().size() + " (" + concat(viruses) + ")";
    }

    @NotNull
    private String wgdString()
    {
        if(PurpleQCInterpretation.isFail(report.purple().fit().qc()))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        return report.purple().characteristics().wholeGenomeDuplication() ? "Yes" : "No";
    }

    @NotNull
    private String msiString()
    {
        if(PurpleQCInterpretation.isFail(report.purple().fit().qc()))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        PurpleCharacteristics characteristics = report.purple().characteristics();
        return formatSingleDigitDecimal(characteristics.microsatelliteIndelsPerMb()) + " ("
                + characteristics.microsatelliteStatus().display() + ")";
    }

    @NotNull
    private String tmbString()
    {
        if(PurpleQCInterpretation.isFail(report.purple().fit().qc()))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        PurpleCharacteristics characteristics = report.purple().characteristics();
        return formatSingleDigitDecimal(characteristics.tumorMutationalBurdenPerMb()) +
                " (" + characteristics.tumorMutationalBurdenStatus().display() + ")";
    }

    @NotNull
    private String tmlString()
    {
        if(PurpleQCInterpretation.isFail(report.purple().fit().qc()))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        PurpleCharacteristics characteristics = report.purple().characteristics();
        return characteristics.tumorMutationalLoad() + " (" + characteristics.tumorMutationalLoadStatus().display() + ")";
    }

    @NotNull
    private String hrDeficiencyString()
    {
        if(PurpleQCInterpretation.isFail(report.purple().fit().qc()))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        ChordRecord chord = report.chord();
        if(chord == null)
        {
            return ReportResources.NOT_AVAILABLE;
        }

        if(chord.hrStatus() == ChordStatus.CANNOT_BE_DETERMINED)
        {
            return chord.hrStatus().display();
        }

        String addon = Strings.EMPTY;
        if(chord.hrStatus() == ChordStatus.HR_DEFICIENT && !chord.hrdType().isEmpty())
        {
            if(chord.hrdType().contains("BRCA1"))
            {
                addon = " - BRCA1 (" + formatTwoDigitDecimal(chord.brca1Value()) + ")";
            }
            else if(chord.hrdType().contains("BRCA2"))
            {
                addon = " - BRCA2 (" + formatTwoDigitDecimal(chord.brca2Value()) + ")";
            }
            else if(chord.hrdType().equals("cannot_be_determined"))
            {
                addon = " - Undetermined";
            }
            else
            {
                addon = " - " + chord.hrdType();
            }
        }
        return formatSingleDigitDecimal(chord.hrdValue()) + " (" + chord.hrStatus().display() + addon + ")";
    }

    @NotNull
    private String geneStatus(@NotNull String gene)
    {
        Set<PeachGenotype> genotypes = report.peach();
        if(genotypes == null)
        {
            return ReportResources.NOT_AVAILABLE;
        }

        Set<String> haplotypes = Sets.newHashSet();
        for(PeachGenotype genotype : genotypes)
        {
            if(genotype.gene().equals(gene))
            {
                String haplotype;
                if(genotype.allele().equals(UNKNOWN_ALLELE_STRING))
                {
                    haplotype = genotype.allele();
                }
                else
                {
                    haplotype = genotype.allele() + " " + convertToZygosityString(genotype.alleleCount());
                }
                haplotypes.add(haplotype + " (" + genotype.function() + ")");
            }
        }
        return !haplotypes.isEmpty() ? concat(haplotypes) : NONE;
    }

    @NotNull
    private String svTmbString()
    {
        if(PurpleQCInterpretation.isFail(report.purple().fit().qc()))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        String svTmb = String.valueOf(report.purple().characteristics().svTumorMutationalBurden());

        Evaluation evaluation = report.cohortEvaluations().get(PercentileType.SV_TMB);
        String addon = Strings.EMPTY;
        if(evaluation != null)
        {
            String panCancerPercentile = formatTwoDigitDecimal(evaluation.panCancerPercentile());
            addon = " (Pan " + panCancerPercentile;
            String cancerType = evaluation.cancerType();
            if(cancerType != null && !cancerType.equals(CohortConstants.COHORT_OTHER)
                    && !cancerType.equals(CohortConstants.COHORT_UNKNOWN))
            {
                Double percentile = evaluation.cancerTypePercentile();
                String cancerTypePercentile = percentile != null ? formatTwoDigitDecimal(percentile) : "NA";
                addon = addon + " | " + evaluation.cancerType() + " " + cancerTypePercentile;
            }
            addon = addon + ")";
        }

        return svTmb + addon;
    }

    @NotNull
    private String maxComplexSizeString()
    {
        if(PurpleQCInterpretation.isFail(report.purple().fit().qc()))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        CuppaData cuppa = report.cuppa();
        return cuppa != null ? Integer.toString(cuppa.maxComplexSize()) : ReportResources.NOT_AVAILABLE;
    }

    @NotNull
    private String telomericSGLString()
    {
        if(PurpleQCInterpretation.isFail(report.purple().fit().qc()))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        CuppaData cuppa = report.cuppa();
        return cuppa != null ? Integer.toString(cuppa.telomericSGLs()) : ReportResources.NOT_AVAILABLE;
    }

    @NotNull
    private String lineCountString()
    {
        if(PurpleQCInterpretation.isFail(report.purple().fit().qc()))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        CuppaData cuppa = report.cuppa();
        return cuppa != null ? Integer.toString(cuppa.lineCount()) : ReportResources.NOT_AVAILABLE;
    }

    @NotNull
    private static String concat(@NotNull Iterable<String> strings)
    {
        StringJoiner joiner = new StringJoiner(", ");
        for(String string : strings)
        {
            joiner.add(string);
        }

        return joiner.toString();
    }
}
