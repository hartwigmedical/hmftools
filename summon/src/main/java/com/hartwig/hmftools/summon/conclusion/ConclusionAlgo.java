package com.hartwig.hmftools.summon.conclusion;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.linx.ReportableHomozygousDisruption;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;
import com.hartwig.hmftools.summon.SummonData;
import com.hartwig.hmftools.summon.actionability.ActionabilityEntry;
import com.hartwig.hmftools.summon.actionability.ActionabilityKey;
import com.hartwig.hmftools.summon.actionability.ImmutableActionabilityKey;
import com.hartwig.hmftools.summon.actionability.Type;

import org.jetbrains.annotations.NotNull;

import com.google.common.collect.Sets;

public class ConclusionAlgo {

    private static final Set<String> FUSION_TYPES = Sets.newHashSet(KnownFusionType.PROMISCUOUS_3.toString(),
            KnownFusionType.PROMISCUOUS_5.toString(),
            KnownFusionType.KNOWN_PAIR.toString(),
            KnownFusionType.IG_KNOWN_PAIR.toString(),
            KnownFusionType.IG_PROMISCUOUS.toString());

    @NotNull
    public static ActionabilityConclusion generateConclusion(@NotNull SummonData summonData) {
        Set<String> conclusion = Sets.newHashSet();

        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = generateActionabilityMap(summonData.actionabilityEntries());

        List<ReportableVariant> reportableSomaticVariants = summonData.purple().reportableSomaticVariants();
        List<ReportableVariant> reportableGermlineVariants = summonData.purple().reportableGermlineVariants();
        List<ReportableGainLoss> reportableGainLosses = summonData.purple().reportableGainsLosses();
        List<LinxFusion> reportableFusions = summonData.linx().reportableFusions();
        List<ReportableHomozygousDisruption> homozygousDisruptions = summonData.linx().homozygousDisruptions();
        List<AnnotatedVirus> reportableViruses = summonData.virusInterpreter().reportableViruses();

        ChordStatus chordStatus = summonData.chord().hrStatus();
        MicrosatelliteStatus microsatelliteStatus = summonData.purple().microsatelliteStatus();
        TumorMutationalStatus tumorMutationalStatus = summonData.purple().tumorMutationalLoadStatus();
        double tumorMutationalBurden = summonData.purple().tumorMutationalBurdenPerMb();

        generateSomaticConclusion(conclusion, reportableSomaticVariants, actionabilityMap);
        generateGermlineConclusion(conclusion, reportableGermlineVariants, actionabilityMap);
        generateCNVConclusion(conclusion, reportableGainLosses, actionabilityMap);
        generateFusionConclusion(conclusion, reportableFusions, actionabilityMap);
        generateHomozygousDisruptionConclusion(conclusion, homozygousDisruptions, actionabilityMap);
        generateVirusConclusion(conclusion, reportableViruses, actionabilityMap);
        generateHrdConclusion(conclusion, chordStatus, actionabilityMap);
        generateMSIConclusion(conclusion, microsatelliteStatus, actionabilityMap);
        generateTMLConclusion(conclusion, tumorMutationalStatus, actionabilityMap);
        generateTMBConclusion(conclusion, tumorMutationalBurden, actionabilityMap);

        return ImmutableActionabilityConclusion.builder().conclusion(conclusion.toString()).build();
    }

    @NotNull
    public static Map<ActionabilityKey, ActionabilityEntry> generateActionabilityMap(@NotNull List<ActionabilityEntry> actionabilityDB) {
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        for (ActionabilityEntry entry : actionabilityDB) {
            ActionabilityKey key = ImmutableActionabilityKey.builder().gene(entry.gene()).type(entry.type()).build();
            actionabilityMap.put(key, entry);
        }
        return actionabilityMap;
    }

    public static void generateSomaticConclusion(@NotNull Set<String> conclusion,
            @NotNull List<ReportableVariant> reportableSomaticVariants,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap) {
        for (ReportableVariant somaticVariant : reportableSomaticVariants) {
            ActionabilityKey keySomaticVariant = ImmutableActionabilityKey.builder()
                    .gene(somaticVariant.gene())
                    .type(somaticVariant.gene().equals("ONCO") ? Type.ACTIVATION_MUTATION : Type.INACTIVATION)
                    .build();
            ActionabilityEntry entry = actionabilityMap.get(keySomaticVariant);
            conclusion.add(entry.conclusion());
        }
    }

    public static void generateGermlineConclusion(@NotNull Set<String> conclusion,
            @NotNull List<ReportableVariant> reportableGermlineVariants,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap) {
        for (ReportableVariant germlineVariant : reportableGermlineVariants) {
            ActionabilityKey keyGermlineVariant = ImmutableActionabilityKey.builder()
                    .gene(germlineVariant.gene())
                    .type(germlineVariant.gene().equals("ONCO") ? Type.ACTIVATION_MUTATION : Type.INACTIVATION)
                    .build();
            ActionabilityEntry entry = actionabilityMap.get(keyGermlineVariant);
            conclusion.add(entry.conclusion());
        }
    }

    public static void generateCNVConclusion(@NotNull Set<String> conclusion, @NotNull List<ReportableGainLoss> reportableGainLosses,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap) {
        for (ReportableGainLoss gainLoss : reportableGainLosses) {
            if (gainLoss.interpretation().display().equals(CopyNumberInterpretation.FULL_LOSS.display())
                    || gainLoss.interpretation().display().equals(CopyNumberInterpretation.PARTIAL_LOSS.display())) {

                ActionabilityKey keyVirus = ImmutableActionabilityKey.builder().gene(gainLoss.gene()).type(Type.LOSS).build();
                ActionabilityEntry entry = actionabilityMap.get(keyVirus);
                conclusion.add(entry.conclusion());
            }

            if (gainLoss.interpretation().display().equals(CopyNumberInterpretation.FULL_GAIN.display())
                    || gainLoss.interpretation().display().equals(CopyNumberInterpretation.PARTIAL_GAIN.display())) {
                ActionabilityKey keyVirus = ImmutableActionabilityKey.builder().gene(gainLoss.gene()).type(Type.AMPLIFICATION).build();
                ActionabilityEntry entry = actionabilityMap.get(keyVirus);
                conclusion.add(entry.conclusion());
            }
        }
    }

    public static void generateFusionConclusion(@NotNull Set<String> conclusion, @NotNull List<LinxFusion> reportableFusions,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap) {
        for (LinxFusion fusion : reportableFusions) {
            if (fusion.reportedType().toString().equals(KnownFusionType.EXON_DEL_DUP.toString())) {
                ActionabilityKey keyFusion =
                        ImmutableActionabilityKey.builder().gene(fusion.geneStart()).type(Type.INTERNAL_DELETION).build();
                ActionabilityEntry entry = actionabilityMap.get(keyFusion);
                conclusion.add(entry.conclusion());
            }
            if (fusion.reportedType().toString().equals(KnownFusionType.EXON_DEL_DUP.toString()) && fusion.geneStart().equals("EGFR")
                    && fusion.fusedExonUp() == 25 && fusion.fusedExonDown() == 18) {
                ActionabilityKey keyFusion =
                        ImmutableActionabilityKey.builder().gene(fusion.geneStart()).type(Type.KINASE_DOMAIN_DUPLICATION).build();
                ActionabilityEntry entry = actionabilityMap.get(keyFusion);
                conclusion.add(entry.conclusion());
            }
            if (FUSION_TYPES.contains(fusion.reportedType().toString())) {
                ActionabilityKey keyFusion = ImmutableActionabilityKey.builder().gene(fusion.geneStart()).type(Type.FUSION).build();
                ActionabilityEntry entry = actionabilityMap.get(keyFusion);
                conclusion.add(entry.conclusion());
            }
        }
    }

    public static void generateHomozygousDisruptionConclusion(@NotNull Set<String> conclusion,
            @NotNull List<ReportableHomozygousDisruption> homozygousDisruptions,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap) {
        for (ReportableHomozygousDisruption homozygousDisruption : homozygousDisruptions) {
            ActionabilityKey keyHomozygousDisruption =
                    ImmutableActionabilityKey.builder().gene(homozygousDisruption.gene()).type(Type.INACTIVATION).build();
            ActionabilityEntry entry = actionabilityMap.get(keyHomozygousDisruption);
            conclusion.add(entry.conclusion());
        }
    }

    public static void generateVirusConclusion(@NotNull Set<String> conclusion, @NotNull List<AnnotatedVirus> reportableViruses,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap) {
        for (AnnotatedVirus annotatedVirus : reportableViruses) {
            if (annotatedVirus.virusDriverLikelihoodType() == VirusLikelihoodType.HIGH) {
                ActionabilityKey keyVirus = ImmutableActionabilityKey.builder().gene("HPV").type(Type.SIGNATURE_POSITIVE).build();
                ActionabilityEntry entry = actionabilityMap.get(keyVirus);
                conclusion.add(entry.conclusion());
            }
        }
    }

    public static void generateHrdConclusion(@NotNull Set<String> conclusion, @NotNull ChordStatus chordStatus,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap) {
        if (chordStatus == ChordStatus.HR_DEFICIENT) {
            ActionabilityKey keyHRD = ImmutableActionabilityKey.builder().gene("HRD").type(Type.SIGNATURE_POSITIVE).build();
            ActionabilityEntry entry = actionabilityMap.get(keyHRD);
            conclusion.add(entry.conclusion());
        }
    }

    public static void generateMSIConclusion(@NotNull Set<String> conclusion, @NotNull MicrosatelliteStatus microsatelliteStatus,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap) {
        if (microsatelliteStatus == MicrosatelliteStatus.MSI) {
            ActionabilityKey keyMSI = ImmutableActionabilityKey.builder().gene("MSI").type(Type.SIGNATURE_POSITIVE).build();
            ActionabilityEntry entry = actionabilityMap.get(keyMSI);
            conclusion.add(entry.conclusion());
        }
    }

    public static void generateTMLConclusion(@NotNull Set<String> conclusion, @NotNull TumorMutationalStatus tumorMutationalStatus,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap) {
        if (tumorMutationalStatus == TumorMutationalStatus.HIGH) {
            ActionabilityKey keyTML = ImmutableActionabilityKey.builder().gene("High-TML").type(Type.SIGNATURE_POSITIVE).build();
            ActionabilityEntry entry = actionabilityMap.get(keyTML);
            conclusion.add(entry.conclusion());
        }
    }

    public static void generateTMBConclusion(@NotNull Set<String> conclusion, double tumorMutationalBurden,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap) {
        if (tumorMutationalBurden >= 10) {
            ActionabilityKey keyTMB = ImmutableActionabilityKey.builder().gene("High-TMB").type(Type.SIGNATURE_POSITIVE).build();
            ActionabilityEntry entry = actionabilityMap.get(keyTMB);
            conclusion.add(entry.conclusion());
        }
    }
}