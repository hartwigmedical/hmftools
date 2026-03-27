package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.peach.PeachUtil.UNKNOWN_ALLELE_STRING;
import static com.hartwig.hmftools.common.peach.PeachUtil.convertToZygosityString;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_NO;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_YES;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatPercentage;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatTwoDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.datamodel.chord.ChordRecord;
import com.hartwig.hmftools.datamodel.chord.ChordStatus;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.cuppa.CuppaPrediction;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.orange.ExperimentType;
import com.hartwig.hmftools.datamodel.orange.OrangeDoidNode;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.datamodel.purple.PurpleCharacteristics;
import com.hartwig.hmftools.datamodel.purple.PurpleFit;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleQC;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.virus.VirusInterpretation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.algo.QcStatusInterpretation;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;

import org.apache.logging.log4j.util.Strings;

public class FrontPageTables
{
    public static Table buildSampleSummary(
            final OrangeRecord report, final OrangeConfig config, float width, final ReportResources reportResources)
    {
        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        // Primary Tumor [as inputed - is it 2 fields?] |  Purity |  Ploidy | Fit Method | QC Status

        String configuredPrimaryTumorLocation = config != null ? config.PrimaryTumorLocation : null;

        String configuredCancerType = configuredPrimaryTumorLocation != null ?
                configuredPrimaryTumorLocation :  configuredPrimaryTumor(report.configuredPrimaryTumor());

        boolean showCancerType = configuredCancerType != null && !configuredCancerType.isEmpty();

        if(showCancerType)
            addEntry(cells, widths, cellEntries, 2, "Primary Tumor");

        addEntry(cells, widths, cellEntries, 2, "Purity");
        addEntry(cells, widths, cellEntries, 2, "Ploidy");
        addEntry(cells, widths, cellEntries, 2, "Fit Method");
        addEntry(cells, widths, cellEntries, 2, "QC");

        Table table = Tables.createContent(width, intToFloatArray(widths), cellArray(cellEntries));

        if(showCancerType)
            table.addCell(cells.createContent(configuredCancerType));

        table.addCell(cells.createContent(purityString(report.purple().fit())));
        table.addCell(cells.createContent(ploidyString(report.purple().fit())));
        table.addCell(cells.createContent(report.purple().fit().fittedPurityMethod().toString()));

        table.addCell(cells.createContent(purpleQCString(report.purple().fit().qc())));

        addQCWarningInCaseOfFail(report.purple().fit().qc(), table, cells);

        return table;
    }

    public static Table buildTechnicalSummary(
            final OrangeRecord report, final OrangeConfig config, float width, final ReportResources reportResources)
    {
        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        // OA Version | Genome Version | Sequencing Type | Pipeline Mode | Panel Name | Date Analysed?

        addEntry(cells, widths, cellEntries, 2, "Pipeline Version");
        addEntry(cells, widths, cellEntries, 2, "Genome Version");
        addEntry(cells, widths, cellEntries, 2, "Sequencing Type");
        addEntry(cells, widths, cellEntries, 2, "Pipeline");
        addEntry(cells, widths, cellEntries, 2, "Samples");
        addEntry(cells, widths, cellEntries, 2, "Date Analysed");

        Table table = Tables.createContent(width, intToFloatArray(widths), cellArray(cellEntries));

        table.addCell(cells.createContent(report.pipelineVersion()));
        table.addCell(cells.createContent(report.refGenomeVersion().toString()));
        table.addCell(cells.createContent(SequencingType.ILLUMINA.toString()));

        // table.addCell(cells.createContent(report.experimentType().toString()));
        table.addCell(cells.createContent(pipelineModeDisplay(report, config)));
        table.addCell(cells.createContent(sampleTypesDisplay(report)));

        table.addCell(cells.createContent(report.samplingDate().toString()));

        return table;
    }

    private static String pipelineModeDisplay(final OrangeRecord report, final OrangeConfig config)
    {
        if(report.experimentType() == ExperimentType.WHOLE_GENOME)
            return "WHOLE GENOME";

        if(config != null && config.PanelName != null)
            return format("PANEL: %s", config.PanelName);
        else
            return format("TARGETED PANEL");
    }

    private static String sampleTypesDisplay(final OrangeRecord report)
    {
        if(report.tumorOnlyMode() && !report.hasRna())
            return "TUMOR-ONLY";

        StringJoiner sj = new StringJoiner(" / ");
        sj.add("TUMOR");

        if(!report.tumorOnlyMode())
            sj.add("NORMAL");

        if(report.hasRna())
            sj.add("RNA");

        return sj.toString();
    }

    private static String configuredPrimaryTumor(final Set<OrangeDoidNode> nodes)
    {
        Set<String> configured = Sets.newHashSet();

        for(OrangeDoidNode node : nodes)
        {
            String primaryTumorInfo = node.doidTerm();

            if(!node.doid().isEmpty())
                primaryTumorInfo += " (DOID " + node.doid() + ")";

            configured.add(primaryTumorInfo);
        }

        return concat(configured);
    }

    private static String purpleQCString(final PurpleQC qc)
    {
        Set<String> purpleStatuses = Sets.newHashSet();
        for(PurpleQCStatus status : qc.status())
        {
            purpleStatuses.add(status.toString());
        }
        return concat(purpleStatuses);
    }

    private static void addQCWarningInCaseOfFail(final PurpleQC qc, final Table table, final Cells cells)
    {
        boolean isFailNoTumor = QcStatusInterpretation.isFailNoTumor(qc);
        boolean isContaminated = QcStatusInterpretation.hasTumorContaminated(qc);

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

    private static String purityString(final PurpleFit purpleFit)
    {
        return format("%s (%s-%s)",
                formatPercentage(purpleFit.purity()), formatPercentage(purpleFit.minPurity()), formatPercentage(purpleFit.maxPurity()));
    }

    private static String ploidyString(final PurpleFit purpleFit)
    {
        return format("%s (%s-%s)",
                formatTwoDigitDecimal(purpleFit.ploidy()), formatTwoDigitDecimal(purpleFit.minPloidy()),
                formatTwoDigitDecimal(purpleFit.maxPloidy()));
    }

    public static Table buildDriverSummary(final OrangeRecord report, float width, final ReportResources reportResources)
    {
        Cells cells = new Cells(reportResources);

        boolean includeGermline = !report.tumorOnlyMode();

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 2 }));

        addCellEntry(table, cells, "Somatic variant:", somaticVariantString(report.purple()));
        addCellEntry(table, cells, "Somatic copy number:", somaticCopyNumberString(report.purple()));
        addCellEntry(table, cells, "Somatic disruption:", somaticDisruptionString(report.linx()));

        if(includeGermline)
        {
            addCellEntry(table, cells, "Germline variant:", germlineVariantString(report.purple()));
            addCellEntry(table, cells, "Germline copy number:", germlineCopyNumberString(report.purple()));
            addCellEntry(table, cells, "Germline disruption:", germlineDisruptionString(report.linx()));
        }

        addCellEntry(table, cells, "Fusion drivers:", fusionsString(report.linx()));

        if(report.virusInterpreter() != null)
        {
            addCellEntry(table, cells, "Viral presence:", virusesString(report.virusInterpreter()));
        }

        addCellEntry(table, cells, "Whole genome duplicated:", wgdString(report.purple()));

        if(includeGermline)
        {
            addCellEntry(table, cells, "DPYD status:", peachGeneStatus(report.peach(), "DPYD"));
            // addCellEntry(table, cells, "UGT1A1 status:", geneStatus("UGT1A1"));
        }

        return new Tables(reportResources).createWrapping(table, "Driver Summary");
    }

    private static final String NONE = "None";

    private static String fusionsString(final LinxRecord linxRecord)
    {
        if(linxRecord.fusions().isEmpty())
            return NONE;

        Set<String> highLikelihoodFusions = Sets.newHashSet();

        for(LinxFusion fusion : linxRecord.fusions())
        {
            if(fusion.driverInterpretation() == DriverInterpretation.LOW)
                continue;

            highLikelihoodFusions.add(fusion.display());
        }

        return reportableEventString(linxRecord.fusions().size(), highLikelihoodFusions);
    }

    private static String reportableEventString(final int driverCount, final Set<String> reportedGenes)
    {
        if(reportedGenes.isEmpty())
            return String.valueOf(driverCount);

        String orderGenesStr = reportedGenes.stream().sorted().collect(Collectors.joining(", "));
        return format("%d (%s)", driverCount, orderGenesStr);
    }

    private static Set<String> findHighDriverVariantGenes(final List<PurpleVariant> variants)
    {
        Set<String> highDriverGenes = Sets.newHashSet();

        for(PurpleVariant variant : variants)
        {
            if(DriverInterpretation.interpret(variant.driverLikelihood()) == DriverInterpretation.LOW)
                continue;

            highDriverGenes.add(variant.gene());
        }

        return highDriverGenes;
    }

    private static String somaticVariantString(final PurpleRecord purpleRecord)
    {
        if(purpleRecord.somaticVariants().isEmpty())
            return NONE;

        Set<String> highDriverGenes = findHighDriverVariantGenes(purpleRecord.somaticVariants());
        return reportableEventString(purpleRecord.somaticVariants().size(), highDriverGenes);
    }

    private static String virusesString(final VirusInterpreterData virusData)
    {
        if(virusData.reportableViruses().isEmpty())
            return NONE;

        Set<String> viruses = Sets.newHashSet();
        for(VirusInterpreterEntry virus : virusData.reportableViruses())
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

        return reportableEventString(virusData.reportableViruses().size(), viruses);
    }

    private static Set<String> findHighDriverAmpDelGenes(final List<PurpleGainDeletion> ampDels)
    {
        Set<String> highDriverGenes = Sets.newHashSet();

        for(PurpleGainDeletion ampDel : ampDels)
        {
            if(DriverInterpretation.interpret(ampDel.driver().driverLikelihood()) == DriverInterpretation.LOW)
                continue;

            highDriverGenes.add(ampDel.gene());
        }

        return highDriverGenes;
    }

    private static String somaticCopyNumberString(final PurpleRecord purpleRecord)
    {
        if(purpleRecord.somaticGainsDels().isEmpty())
            return NONE;

        Set<String> highDriverGenes = findHighDriverAmpDelGenes(purpleRecord.somaticGainsDels());
        return reportableEventString(purpleRecord.somaticGainsDels().size(), highDriverGenes);
    }

    private static Set<String> findHighDriverDisruptionGenes(final List<LinxBreakend> breakends)
    {
        Set<String> highDriverGenes = Sets.newHashSet();

        for(LinxBreakend breakend : breakends)
        {
            if(DriverInterpretation.interpret(breakend.driverLikelihood()) == DriverInterpretation.LOW)
                continue;

            highDriverGenes.add(breakend.gene());
        }

        return highDriverGenes;
    }

    private static String somaticDisruptionString(final LinxRecord linxRecord)
    {
        if(linxRecord.somaticBreakends().isEmpty())
            return NONE;

        Set<String> highDriverGenes = findHighDriverDisruptionGenes(linxRecord.somaticBreakends());
        return reportableEventString(linxRecord.somaticBreakends().size(), highDriverGenes);
    }

    private static String germlineVariantString(final PurpleRecord purpleRecord)
    {
        if(purpleRecord.germlineVariants() == null || purpleRecord.germlineVariants().isEmpty())
            return NONE;

        Set<String> highDriverGenes = findHighDriverVariantGenes(purpleRecord.germlineVariants());
        return reportableEventString(purpleRecord.germlineVariants().size(), highDriverGenes);

    }

    private static String germlineCopyNumberString(final PurpleRecord purpleRecord)
    {
        if(purpleRecord.germlineGainsDels() == null || purpleRecord.germlineGainsDels().isEmpty())
            return NONE;

        Set<String> highDriverGenes = findHighDriverAmpDelGenes(purpleRecord.germlineGainsDels());
        return reportableEventString(purpleRecord.germlineGainsDels().size(), highDriverGenes);
    }

    private static String germlineDisruptionString(final LinxRecord linxRecord)
    {
        if(linxRecord.germlineBreakends() == null || linxRecord.germlineBreakends().isEmpty())
            return NONE;

        Set<String> highDriverGenes = findHighDriverDisruptionGenes(linxRecord.germlineBreakends());
        return reportableEventString(linxRecord.germlineBreakends().size(), highDriverGenes);
    }

    public static Table buildGenomeWideFeatures(final OrangeRecord report, float width, final ReportResources reportResources)
    {
        boolean includeGermline = !report.tumorOnlyMode();

        Cells cells = new Cells(reportResources);

        float labelWidth = (float)(1.9);
        float valueWidth = (float)(1.0);
        Table table = new Table(UnitValue.createPercentArray(new float[] { labelWidth, valueWidth }));

        addCellEntry(table, cells, "Microsatellite indels per Mb:", msiString(report.purple()));
        addCellEntry(table, cells, "Tumor mutations per Mb:", tmbString(report.purple()));
        addCellEntry(table, cells, "Tumor mutational load:", tmlString(report.purple()));

        if(includeGermline) // will change once we solve HRD for targeted mode
        {
            addCellEntry(table, cells, "HR deficiency score:", hrDeficiencyString(report.chord()));
        }

        addCellEntry(table, cells, "LOH proportion:", lohPercentageString(report.purple()));
        addCellEntry(table, cells, "Number of SVs:", svTmbString(report.purple()));

        if(!report.tumorOnlyMode() && report.cuppa() != null)
            addCellEntry(table, cells, "CUPPA cancer type:", cuppaCancerType(report.cuppa()));

        return new Tables(reportResources).createWrapping(table, "Genome Wide Biomarkers");
    }

    private static String cuppaCancerType(final CuppaData cuppaData)
    {
        CuppaPrediction best = cuppaData.bestPrediction();
        return best.cancerType() + " (" + formatPercentage(best.likelihood()) + ")";
    }

    private static String wgdString(final PurpleRecord purpleRecord)
    {
        return purpleRecord.characteristics().wholeGenomeDuplication() ? VALUE_YES : VALUE_NO;
    }

    private static String msiString(final PurpleRecord purpleRecord)
    {
        PurpleCharacteristics characteristics = purpleRecord.characteristics();
        return formatSingleDigitDecimal(characteristics.microsatelliteIndelsPerMb()) + " ("
                + characteristics.microsatelliteStatus().display() + ")";
    }

    private static String tmbString(final PurpleRecord purpleRecord)
    {
        PurpleCharacteristics characteristics = purpleRecord.characteristics();
        return formatSingleDigitDecimal(characteristics.tumorMutationalBurdenPerMb()) +
                " (" + characteristics.tumorMutationalBurdenStatus().display() + ")";
    }

    private static String tmlString(final PurpleRecord purpleRecord)
    {
        PurpleCharacteristics characteristics = purpleRecord.characteristics();
        return characteristics.tumorMutationalLoad() + " (" + characteristics.tumorMutationalLoadStatus().display() + ")";
    }

    private static String hrDeficiencyString(final ChordRecord chord)
    {
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

    private static String peachGeneStatus(final Set<PeachGenotype> peachGenotypes, final String gene)
    {
        if(peachGenotypes == null)
        {
            return ReportResources.NOT_AVAILABLE;
        }

        Set<String> haplotypes = Sets.newHashSet();
        for(PeachGenotype genotype : peachGenotypes)
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
        return !haplotypes.isEmpty() ? concat(haplotypes) : ReportResources.NOT_AVAILABLE;
    }

    private static String svTmbString(final PurpleRecord purpleRecord)
    {
        String svTmb = String.valueOf(purpleRecord.characteristics().svTumorMutationalBurden());
        String addon = Strings.EMPTY;
        return svTmb + addon;
    }

    private static String lohPercentageString(final PurpleRecord purpleRecord)
    {
        String lohPerc = formatPercentage(purpleRecord.characteristics().lohPercentage());
        String addon = Strings.EMPTY;
        return lohPerc + addon;
    }

    private static void addCellEntry(final Table summary, final Cells cells, final String name, final String value)
    {
        summary.addCell(cells.createKey(name));
        summary.addCell(cells.createValue(value));
    }

    private static String concat(final Iterable<String> strings)
    {
        StringJoiner joiner = new StringJoiner(", ");
        for(String string : strings)
        {
            joiner.add(string);
        }

        return joiner.toString();
    }
}
