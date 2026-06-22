package com.hartwig.hmftools.orange.report.pdfdata;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.peach.PeachUtil.UNKNOWN_ALLELE_STRING;
import static com.hartwig.hmftools.common.peach.PeachUtil.convertToZygosityString;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatPercentage;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatTwoDigitDecimal;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

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
import com.hartwig.hmftools.orange.report.PlotPathResolver;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.Nullable;

public class FrontPageDataFactory
{
    private static final String NONE = "None";

    public static FrontPageData build(final OrangeRecord report, @Nullable final OrangeConfig config,
            final PlotPathResolver plotPathResolver)
    {
        return new FrontPageData(
                buildSampleSummary(report, config),
                buildQcWarning(report),
                buildTechnicalSummary(report, config),
                buildDriverSummary(report),
                buildGenomeWideFeatures(report),
                plotPathResolver.resolve(report.plots().purpleFinalCircosPlot()));
    }

    private static Map<String, String> buildSampleSummary(final OrangeRecord report, @Nullable final OrangeConfig config)
    {
        String configuredPrimaryTumorLocation = config != null ? config.PrimaryTumorLocation : null;
        String configuredCancerType = configuredPrimaryTumorLocation != null
                ? configuredPrimaryTumorLocation
                : configuredPrimaryTumor(report.configuredPrimaryTumor());
        boolean hasCancerType = configuredCancerType != null && !configuredCancerType.isEmpty();

        Map<String, String> summary = new LinkedHashMap<>();
        summary.put("Primary Tumor", hasCancerType ? configuredCancerType : "NOT SPECIFIED");
        summary.put("Purity", purityString(report.purple().fit()));
        summary.put("Ploidy", ploidyString(report.purple().fit()));
        summary.put("Fit Method", report.purple().fit().fittedPurityMethod().toString());
        summary.put("QC", purpleQCString(report.purple().fit().qc()));
        return summary;
    }

    @Nullable
    private static String buildQcWarning(final OrangeRecord report)
    {
        PurpleQC qc = report.purple().fit().qc();
        boolean isFailNoTumor = QcStatusInterpretation.isFailNoTumor(qc);
        boolean isContaminated = QcStatusInterpretation.hasTumorContaminated(qc);

        if(!isFailNoTumor && !isContaminated)
        {
            return null;
        }

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

        return "The QC status of this sample is fail (" + reason + ")"
                + ": all presented data in this report should be interpreted with caution";
    }

    private static Map<String, String> buildTechnicalSummary(final OrangeRecord report, @Nullable final OrangeConfig config)
    {
        Map<String, String> summary = new LinkedHashMap<>();
        summary.put("Pipeline Version", report.pipelineVersion());
        summary.put("Genome Version", report.refGenomeVersion().toString());
        summary.put("Sequencing Type", sequencingType(config));
        summary.put("Pipeline", pipelineModeDisplay(report, config));
        summary.put("Samples", sampleTypesDisplay(report));
        summary.put("Date Analysed", report.samplingDate().toString());
        return summary;
    }

    private static Map<String, String> buildDriverSummary(final OrangeRecord report)
    {
        boolean includeGermline = !report.tumorOnlyMode();

        Map<String, String> summary = new LinkedHashMap<>();
        summary.put("Somatic variant:", somaticVariantString(report.purple()));
        summary.put("Somatic copy number:", somaticCopyNumberString(report.purple()));
        summary.put("Somatic disruption:", somaticDisruptionString(report.linx()));

        if(includeGermline)
        {
            summary.put("Germline variant:", germlineVariantString(report.purple()));
            summary.put("Germline copy number:", germlineCopyNumberString(report.purple()));
            summary.put("Germline disruption:", germlineDisruptionString(report.linx()));
        }

        summary.put("Fusion drivers:", fusionsString(report.linx()));

        if(report.virusInterpreter() != null)
        {
            summary.put("Viral presence:", virusesString(report.virusInterpreter()));
        }

        summary.put("Whole genome duplicated:", report.purple().characteristics().wholeGenomeDuplication() ? "Yes" : "No");

        if(includeGermline)
        {
            summary.put("DPYD status:", peachGeneStatus(report.peach(), "DPYD"));
        }

        return summary;
    }

    private static Map<String, String> buildGenomeWideFeatures(final OrangeRecord report)
    {
        boolean includeGermline = !report.tumorOnlyMode();

        Map<String, String> features = new LinkedHashMap<>();
        features.put("Microsatellite indels per Mb:", msiString(report.purple()));
        features.put("Tumor mutations per Mb:", tmbString(report.purple()));
        features.put("Tumor mutational load:", tmlString(report.purple()));

        if(includeGermline)
        {
            features.put("HR deficiency score:", hrDeficiencyString(report.chord()));
        }

        features.put("LOH proportion:", lohPercentageString(report.purple()));
        features.put("Number of SVs:", svTmbString(report.purple()));

        if(!report.tumorOnlyMode() && report.cuppa() != null)
        {
            features.put("CUPPA cancer type:", cuppaCancerType(report.cuppa()));
        }

        return features;
    }

    private static String configuredPrimaryTumor(final Set<OrangeDoidNode> nodes)
    {
        Set<String> configured = Sets.newHashSet();
        for(OrangeDoidNode node : nodes)
        {
            String primaryTumorInfo = node.doidTerm();
            if(!node.doid().isEmpty())
            {
                primaryTumorInfo += " (DOID " + node.doid() + ")";
            }
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

    private static String sequencingType(@Nullable final OrangeConfig config)
    {
        return config != null ? config.SeqType.toString() : SequencingType.ILLUMINA.toString();
    }

    private static String pipelineModeDisplay(final OrangeRecord report, @Nullable final OrangeConfig config)
    {
        if(report.experimentType() == ExperimentType.WHOLE_GENOME)
        {
            return "WHOLE GENOME";
        }

        if(config != null && config.PanelName != null)
        {
            return format("PANEL: %s", config.PanelName);
        }
        else
        {
            return "TARGETED PANEL";
        }
    }

    private static String sampleTypesDisplay(final OrangeRecord report)
    {
        if(report.tumorOnlyMode() && !report.hasRna())
        {
            return "TUMOR-ONLY";
        }

        StringJoiner sj = new StringJoiner(" / ");
        sj.add("TUMOR");
        if(!report.tumorOnlyMode())
        {
            sj.add("NORMAL");
        }
        if(report.hasRna())
        {
            sj.add("RNA");
        }
        return sj.toString();
    }

    private static String fusionsString(final LinxRecord linxRecord)
    {
        if(linxRecord.fusions().isEmpty())
        {
            return NONE;
        }

        Set<String> highLikelihoodFusions = Sets.newHashSet();
        for(LinxFusion fusion : linxRecord.fusions())
        {
            if(fusion.driverInterpretation() == DriverInterpretation.LOW)
            {
                continue;
            }
            highLikelihoodFusions.add(fusion.display());
        }

        return reportableEventString(linxRecord.fusions().size(), highLikelihoodFusions);
    }

    private static String reportableEventString(final int driverCount, final Set<String> reportedGenes)
    {
        if(reportedGenes.isEmpty())
        {
            return String.valueOf(driverCount);
        }

        String orderedGenesStr = reportedGenes.stream().sorted().collect(Collectors.joining(", "));
        return format("%d (%s)", driverCount, orderedGenesStr);
    }

    private static Set<String> findHighDriverVariantGenes(final List<PurpleVariant> variants)
    {
        Set<String> highDriverGenes = Sets.newHashSet();
        for(PurpleVariant variant : variants)
        {
            if(DriverInterpretation.interpret(variant.driverLikelihood()) == DriverInterpretation.LOW)
            {
                continue;
            }
            highDriverGenes.add(variant.gene());
        }
        return highDriverGenes;
    }

    private static String somaticVariantString(final PurpleRecord purpleRecord)
    {
        if(purpleRecord.somaticVariants().isEmpty())
        {
            return NONE;
        }

        Set<String> highDriverGenes = findHighDriverVariantGenes(purpleRecord.somaticVariants());
        return reportableEventString(purpleRecord.somaticVariants().size(), highDriverGenes);
    }

    private static String virusesString(final VirusInterpreterData virusData)
    {
        if(virusData.reportableViruses().isEmpty())
        {
            return NONE;
        }

        Set<String> viruses = Sets.newHashSet();
        for(VirusInterpreterEntry virus : virusData.reportableViruses())
        {
            VirusInterpretation interpretation = virus.interpretation();
            viruses.add(interpretation != null ? interpretation.name() : virus.name());
        }

        return reportableEventString(virusData.reportableViruses().size(), viruses);
    }

    private static Set<String> findHighDriverAmpDelGenes(final List<PurpleGainDeletion> ampDels)
    {
        Set<String> highDriverGenes = Sets.newHashSet();
        for(PurpleGainDeletion ampDel : ampDels)
        {
            if(DriverInterpretation.interpret(ampDel.driver().driverLikelihood()) == DriverInterpretation.LOW)
            {
                continue;
            }
            highDriverGenes.add(ampDel.gene());
        }
        return highDriverGenes;
    }

    private static String somaticCopyNumberString(final PurpleRecord purpleRecord)
    {
        if(purpleRecord.somaticGainsDels().isEmpty())
        {
            return NONE;
        }

        Set<String> highDriverGenes = findHighDriverAmpDelGenes(purpleRecord.somaticGainsDels());
        return reportableEventString(purpleRecord.somaticGainsDels().size(), highDriverGenes);
    }

    private static Set<String> findHighDriverDisruptionGenes(final List<LinxBreakend> breakends)
    {
        Set<String> highDriverGenes = Sets.newHashSet();
        for(LinxBreakend breakend : breakends)
        {
            if(DriverInterpretation.interpret(breakend.driverLikelihood()) == DriverInterpretation.LOW)
            {
                continue;
            }
            highDriverGenes.add(breakend.gene());
        }
        return highDriverGenes;
    }

    private static String somaticDisruptionString(final LinxRecord linxRecord)
    {
        if(linxRecord.somaticBreakends().isEmpty())
        {
            return NONE;
        }

        Set<String> highDriverGenes = findHighDriverDisruptionGenes(linxRecord.somaticBreakends());
        return reportableEventString(linxRecord.somaticBreakends().size(), highDriverGenes);
    }

    private static String germlineVariantString(final PurpleRecord purpleRecord)
    {
        if(purpleRecord.germlineVariants() == null || purpleRecord.germlineVariants().isEmpty())
        {
            return NONE;
        }

        Set<String> highDriverGenes = findHighDriverVariantGenes(purpleRecord.germlineVariants());
        return reportableEventString(purpleRecord.germlineVariants().size(), highDriverGenes);
    }

    private static String germlineCopyNumberString(final PurpleRecord purpleRecord)
    {
        if(purpleRecord.germlineGainsDels() == null || purpleRecord.germlineGainsDels().isEmpty())
        {
            return NONE;
        }

        Set<String> highDriverGenes = findHighDriverAmpDelGenes(purpleRecord.germlineGainsDels());
        return reportableEventString(purpleRecord.germlineGainsDels().size(), highDriverGenes);
    }

    private static String germlineDisruptionString(final LinxRecord linxRecord)
    {
        if(linxRecord.germlineBreakends() == null || linxRecord.germlineBreakends().isEmpty())
        {
            return NONE;
        }

        Set<String> highDriverGenes = findHighDriverDisruptionGenes(linxRecord.germlineBreakends());
        return reportableEventString(linxRecord.germlineBreakends().size(), highDriverGenes);
    }

    private static String cuppaCancerType(final CuppaData cuppaData)
    {
        CuppaPrediction best = cuppaData.bestPrediction();
        return best.cancerType() + " (" + formatPercentage(best.likelihood()) + ")";
    }

    private static String msiString(final PurpleRecord purpleRecord)
    {
        PurpleCharacteristics characteristics = purpleRecord.characteristics();
        return formatSingleDigitDecimal(characteristics.microsatelliteIndelsPerMb())
                + " (" + characteristics.microsatelliteStatus().display() + ")";
    }

    private static String tmbString(final PurpleRecord purpleRecord)
    {
        PurpleCharacteristics characteristics = purpleRecord.characteristics();
        return formatSingleDigitDecimal(characteristics.tumorMutationalBurdenPerMb())
                + " (" + characteristics.tumorMutationalBurdenStatus().display() + ")";
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
            if(!genotype.gene().equals(gene))
            {
                continue;
            }

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
        return !haplotypes.isEmpty() ? concat(haplotypes) : ReportResources.NOT_AVAILABLE;
    }

    private static String svTmbString(final PurpleRecord purpleRecord)
    {
        return String.valueOf(purpleRecord.characteristics().svTumorMutationalBurden());
    }

    private static String lohPercentageString(final PurpleRecord purpleRecord)
    {
        return formatPercentage(purpleRecord.characteristics().lohPercentage());
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
