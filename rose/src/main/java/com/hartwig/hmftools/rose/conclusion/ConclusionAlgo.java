package com.hartwig.hmftools.rose.conclusion;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumorFunctions;
import com.hartwig.hmftools.common.cuppa.MolecularTissueOrigin;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.linx.HomozygousDisruption;
import com.hartwig.hmftools.common.protect.EventGenerator;
import com.hartwig.hmftools.common.purple.loader.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.loader.GainLoss;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.DriverInterpretation;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantFactory;
import com.hartwig.hmftools.common.variant.ReportableVariantSource;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.purple.TumorMutationalStatus;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.rose.RoseData;
import com.hartwig.hmftools.rose.actionability.ActionabilityEntry;
import com.hartwig.hmftools.rose.actionability.ActionabilityKey;
import com.hartwig.hmftools.rose.actionability.Condition;
import com.hartwig.hmftools.rose.actionability.ImmutableActionabilityKey;
import com.hartwig.hmftools.rose.actionability.TypeAlteration;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ConclusionAlgo {

    private static final Set<String> FUSION_TYPES = Sets.newHashSet(KnownFusionType.PROMISCUOUS_3.toString(),
            KnownFusionType.PROMISCUOUS_5.toString(),
            KnownFusionType.KNOWN_PAIR.toString(),
            KnownFusionType.IG_KNOWN_PAIR.toString(),
            KnownFusionType.IG_PROMISCUOUS.toString());
    private static final Set<String> HRD_GENES = Sets.newHashSet("BRCA1", "BRCA2", "PALB2", "RAD51B", "RAD51C");

    private static final DecimalFormat DOUBLE_DECIMAL_FORMAT = decimalFormat("#.##");
    private static final double TMB_CUTOFF = 10;
    private static final double PURITY_CUTOFF = 0.195;

    @NotNull
    public static ActionabilityConclusion generateConclusion(@NotNull RoseData roseData) {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Set<String> oncogenic = Sets.newHashSet();
        Set<String> actionable = Sets.newHashSet();
        Set<String> HRD = Sets.newHashSet();

        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = generateActionabilityMap(roseData.actionabilityEntries());
        Map<String, DriverGene> driverGenesMap = generateDriverGenesMap(roseData.driverGenes());
        List<ReportableVariant> reportableSomaticVariants = roseData.purple().reportableSomaticVariants();
        List<ReportableVariant> reportableGermlineVariants = roseData.purple().reportableGermlineVariants();
        List<ReportableVariant> reportableVariants =
                ReportableVariantFactory.mergeVariantLists(reportableGermlineVariants, reportableSomaticVariants);
        List<GainLoss> reportableGainLosses = roseData.purple().reportableSomaticGainsLosses();
        List<LinxFusion> reportableFusions = roseData.linx().reportableFusions();
        List<HomozygousDisruption> homozygousDisruptions = roseData.linx().homozygousDisruptions();
        List<AnnotatedVirus> reportableViruses = roseData.virusInterpreter().reportableViruses();

        genertateTumorLocationConclusion(conclusion, roseData.patientPrimaryTumors(), roseData.patientId());
        genertatePurityConclusion(conclusion, roseData.purple().purity(), roseData.purple().hasReliablePurity(), actionabilityMap);
        generateCUPPAConclusion(conclusion, roseData.molecularTissueOrigin(), actionabilityMap);
        generateVariantConclusion(conclusion, reportableVariants, actionabilityMap, driverGenesMap, oncogenic, actionable, HRD);
        generateCNVConclusion(conclusion, reportableGainLosses, actionabilityMap, oncogenic, actionable);
        generateFusionConclusion(conclusion, reportableFusions, actionabilityMap, oncogenic, actionable);
        generateHomozygousDisruptionConclusion(conclusion, homozygousDisruptions, actionabilityMap, oncogenic, actionable);
        generateVirusConclusion(conclusion, reportableViruses, actionabilityMap, oncogenic, actionable);
        generateHrdConclusion(conclusion, roseData.chord(), actionabilityMap, oncogenic, actionable, HRD);
        generateMSIConclusion(conclusion,
                roseData.purple().microsatelliteStatus(),
                roseData.purple().microsatelliteIndelsPerMb(),
                actionabilityMap,
                oncogenic,
                actionable);
        generateTMLConclusion(conclusion,
                roseData.purple().tumorMutationalLoadStatus(),
                roseData.purple().tumorMutationalLoad(),
                actionabilityMap,
                oncogenic,
                actionable);
        generateTMBConclusion(conclusion, roseData.purple().tumorMutationalBurdenPerMb(), actionabilityMap, oncogenic, actionable);

        generateTotalResults(conclusion, actionabilityMap, oncogenic, actionable);
        generateFindings(conclusion, actionabilityMap);

        String conclusionString = generateConslusionString(conclusion);

        return ImmutableActionabilityConclusion.builder().conclusion(conclusionString).build();
    }

    @NotNull
    public static String generateConslusionString(@NotNull Map<Integer, String> conclusion) {
        StringBuilder conclusionBuilder = new StringBuilder();
        int location = 0;
        for (Map.Entry<Integer, String> entry : conclusion.entrySet()) {
            if (entry.getKey().equals(location)) {
                conclusionBuilder.append(entry.getValue()).append(" <enter> ");
                location += 1;
            }
        }
        return conclusionBuilder.toString();
    }

    @NotNull
    public static Map<ActionabilityKey, ActionabilityEntry> generateActionabilityMap(@NotNull List<ActionabilityEntry> actionabilityDB) {
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        for (ActionabilityEntry entry : actionabilityDB) {
            ActionabilityKey key = ImmutableActionabilityKey.builder().match(entry.match()).type(entry.type()).build();
            actionabilityMap.put(key, entry);
        }
        return actionabilityMap;
    }

    @NotNull
    public static Map<String, DriverGene> generateDriverGenesMap(@NotNull List<DriverGene> driverGenes) {
        Map<String, DriverGene> driverGeneMap = Maps.newHashMap();
        for (DriverGene entry : driverGenes) {
            driverGeneMap.put(entry.gene(), entry);
        }
        return driverGeneMap;
    }

    public static void genertateTumorLocationConclusion(@NotNull Map<Integer, String> conclusion,
            @NotNull List<PatientPrimaryTumor> patientPrimaryTumors, @NotNull String patientId) {
        PatientPrimaryTumor patientPrimaryTumor = PatientPrimaryTumorFunctions.findPrimaryTumorForPatient(patientPrimaryTumors, patientId);

        conclusion.put(conclusion.size(), resolveTumorLocation(patientPrimaryTumor) + " sample showing: ");
    }

    @NotNull
    @VisibleForTesting
    static String resolveTumorLocation(@Nullable PatientPrimaryTumor patientPrimaryTumor) {
        String mergedLocation = Strings.EMPTY;
        String mergedType = Strings.EMPTY;

        if (patientPrimaryTumor != null) {
            if (!patientPrimaryTumor.subLocation().isEmpty()) {
                mergedLocation = patientPrimaryTumor.subLocation();
            } else {
                mergedLocation = patientPrimaryTumor.location();
            }

            if (!patientPrimaryTumor.subType().isEmpty()) {
                mergedType = patientPrimaryTumor.subType();
            } else {
                mergedType = patientPrimaryTumor.type();
            }
        }

        if (!mergedType.isEmpty()) {
            return mergedLocation + " (" + mergedType + ")";
        } else {
            return mergedLocation;
        }
    }

    public static void generateCUPPAConclusion(@NotNull Map<Integer, String> conclusion, MolecularTissueOrigin molecularTissueOrigin,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap) {

        if (molecularTissueOrigin.conclusion().contains("results inconclusive")) {
            ActionabilityKey keyCuppaInconclusice =
                    ImmutableActionabilityKey.builder().match("CUPPA_INCONCLUSIVE").type(TypeAlteration.CUPPA_INCONCLUSIVE).build();

            ActionabilityEntry entry = actionabilityMap.get(keyCuppaInconclusice);
            if (entry != null && entry.condition() == Condition.OTHER) {
                conclusion.put(conclusion.size(), "- " + entry.conclusion());
            }
        } else {
            ActionabilityKey keyCuppa = ImmutableActionabilityKey.builder().match("CUPPA").type(TypeAlteration.CUPPA).build();

            ActionabilityEntry entry = actionabilityMap.get(keyCuppa);
            if (entry != null && entry.condition() == Condition.OTHER) {
                conclusion.put(conclusion.size(), "- " + entry.conclusion() + " " + molecularTissueOrigin.conclusion());
            }
        }
    }

    public static void generateVariantConclusion(@NotNull Map<Integer, String> conclusion,
            @NotNull List<ReportableVariant> reportableVariants, @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap,
            @NotNull Map<String, DriverGene> driverGenesMap, @NotNull Set<String> oncogenic, @NotNull Set<String> actionable,
            @NotNull Set<String> HRD) {

        for (ReportableVariant reportableVariant : reportableVariants) {
            String variant = EventGenerator.variantEvent(reportableVariant);

            if (HRD_GENES.contains(reportableVariant.gene())) {
                HRD.add(reportableVariant.gene());
            }
            oncogenic.add(reportableVariant.source() == ReportableVariantSource.SOMATIC ? "somaticVariant" : "germlineVariant");
            ActionabilityKey keySomaticVariant = ImmutableActionabilityKey.builder()
                    .match(reportableVariant.gene())
                    .type(driverGenesMap.get(reportableVariant.gene()).likelihoodType().equals(DriverCategory.ONCO)
                            ? TypeAlteration.ACTIVATING_MUTATION
                            : TypeAlteration.INACTIVATION)
                    .build();
            ActionabilityEntry entry = actionabilityMap.get(keySomaticVariant);
            if (entry != null) {
                if ((reportableVariant.driverLikelihoodInterpretation() == DriverInterpretation.HIGH
                        && entry.condition() == Condition.ONLY_HIGH) || entry.condition() == Condition.ALWAYS_NO_ACTIONABLE) {
                    if (entry.condition() == Condition.ONLY_HIGH) {
                        actionable.add(
                                reportableVariant.source() == ReportableVariantSource.SOMATIC ? "somaticVariant" : "germlineVariant");
                    }
                    //TODO: Add sentence about germline findings probably in future
                    if (driverGenesMap.get(reportableVariant.gene()).likelihoodType().equals(DriverCategory.TSG)
                            && !reportableVariant.biallelic()) {
                        ActionabilityKey keyBiallelic =
                                ImmutableActionabilityKey.builder().match("NOT_BIALLELIC").type(TypeAlteration.NOT_BIALLELIC).build();
                        ActionabilityEntry entryBiallelic = actionabilityMap.get(keyBiallelic);
                        if (entryBiallelic.condition() == Condition.OTHER) {
                            conclusion.put(conclusion.size(),
                                    "- " + reportableVariant.gene() + " (" + variant + ") " + entry.conclusion() + " "
                                            + entryBiallelic.conclusion());
                        }
                    } else {
                        conclusion.put(conclusion.size(), "- " + reportableVariant.gene() + " (" + variant + ") " + entry.conclusion());
                    }
                }
            }
        }
    }

    public static void generateCNVConclusion(@NotNull Map<Integer, String> conclusion, @NotNull List<GainLoss> reportableGainLosses,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap, @NotNull Set<String> oncogenic,
            @NotNull Set<String> actionable) {
        for (GainLoss gainLoss : reportableGainLosses) {
            oncogenic.add("CNV");

            if (gainLoss.interpretation().display().equals(CopyNumberInterpretation.FULL_LOSS.display()) || gainLoss.interpretation()
                    .display()
                    .equals(CopyNumberInterpretation.PARTIAL_LOSS.display())) {

                ActionabilityKey keyVirus = ImmutableActionabilityKey.builder().match(gainLoss.gene()).type(TypeAlteration.LOSS).build();
                ActionabilityEntry entry = actionabilityMap.get(keyVirus);

                if (entry != null && (entry.condition() == Condition.ALWAYS || entry.condition() == Condition.ALWAYS_NO_ACTIONABLE)) {
                    conclusion.put(conclusion.size(), "- " + gainLoss.gene() + " " + entry.conclusion());
                    actionable.add("CNV");
                }
            }

            if (gainLoss.interpretation().display().equals(CopyNumberInterpretation.FULL_GAIN.display()) || gainLoss.interpretation()
                    .display()
                    .equals(CopyNumberInterpretation.PARTIAL_GAIN.display())) {
                ActionabilityKey keyVirus =
                        ImmutableActionabilityKey.builder().match(gainLoss.gene()).type(TypeAlteration.AMPLIFICATION).build();
                ActionabilityEntry entry = actionabilityMap.get(keyVirus);

                if (entry != null && entry.condition() == Condition.ALWAYS) {
                    conclusion.put(conclusion.size(), "- " + gainLoss.gene() + " " + entry.conclusion());
                    actionable.add("CNV");
                }
            }
        }
    }

    public static void generateFusionConclusion(@NotNull Map<Integer, String> conclusion, @NotNull List<LinxFusion> reportableFusions,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap, @NotNull Set<String> oncogenic,
            @NotNull Set<String> actionable) {
        for (LinxFusion fusion : reportableFusions) {
            oncogenic.add("fusion");

            if (fusion.reportedType().equals(KnownFusionType.EXON_DEL_DUP.toString()) && fusion.geneStart().equals("EGFR") && (
                    fusion.fusedExonUp() == 25 && fusion.fusedExonDown() == 14) || (fusion.fusedExonUp() == 26
                    && fusion.fusedExonDown() == 18)) {
                ActionabilityKey keyFusion = ImmutableActionabilityKey.builder()
                        .match(fusion.geneStart())
                        .type(TypeAlteration.KINASE_DOMAIN_DUPLICATION)
                        .build();
                ActionabilityEntry entry = actionabilityMap.get(keyFusion);
                if (entry != null && entry.condition() == Condition.ALWAYS) {
                    conclusion.put(conclusion.size(), "- " + fusion.name() + " " + entry.conclusion());
                    actionable.add("fusion");
                }
            } else if (fusion.reportedType().equals(KnownFusionType.EXON_DEL_DUP.toString())) {
                ActionabilityKey keyFusion =
                        ImmutableActionabilityKey.builder().match(fusion.geneStart()).type(TypeAlteration.INTERNAL_DELETION).build();
                ActionabilityEntry entry = actionabilityMap.get(keyFusion);
                if (entry != null && entry.condition() == Condition.ALWAYS) {
                    conclusion.put(conclusion.size(), "- " + fusion.name() + " " + entry.conclusion());
                    actionable.add("fusion");
                }
            } else if (FUSION_TYPES.contains(fusion.reportedType())) {
                ActionabilityKey keyFusionStart =
                        ImmutableActionabilityKey.builder().match(fusion.geneStart()).type(TypeAlteration.FUSION).build();
                ActionabilityKey keyFusionEnd =
                        ImmutableActionabilityKey.builder().match(fusion.geneEnd()).type(TypeAlteration.FUSION).build();

                ActionabilityEntry entryStart = actionabilityMap.get(keyFusionStart);
                ActionabilityEntry entryEnd = actionabilityMap.get(keyFusionEnd);

                if (entryStart != null && entryStart.condition() == Condition.ALWAYS) {
                    conclusion.put(conclusion.size(), "- " + fusion.name() + " " + entryStart.conclusion());
                    actionable.add("fusion");
                } else if (entryEnd != null && entryEnd.condition() == Condition.ALWAYS) {
                    conclusion.put(conclusion.size(), "- " + fusion.name() + " " + entryEnd.conclusion());
                    actionable.add("fusion");
                }
            }
        }
    }

    public static void generateHomozygousDisruptionConclusion(@NotNull Map<Integer, String> conclusion,
            @NotNull List<HomozygousDisruption> homozygousDisruptions,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap, @NotNull Set<String> oncogenic,
            @NotNull Set<String> actionable) {
        for (HomozygousDisruption homozygousDisruption : homozygousDisruptions) {
            oncogenic.add("homozygousDisruption");

            ActionabilityKey keyHomozygousDisruption =
                    ImmutableActionabilityKey.builder().match(homozygousDisruption.gene()).type(TypeAlteration.INACTIVATION).build();
            ActionabilityEntry entry = actionabilityMap.get(keyHomozygousDisruption);
            if (entry != null && entry.condition() == Condition.ALWAYS) {
                conclusion.put(conclusion.size(), "- " + homozygousDisruption.gene() + " " + entry.conclusion());
                actionable.add("homozygousDisruption");
            }
        }
    }

    public static void generateVirusConclusion(@NotNull Map<Integer, String> conclusion, @NotNull List<AnnotatedVirus> reportableViruses,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap, @NotNull Set<String> oncogenic,
            @NotNull Set<String> actionable) {
        for (AnnotatedVirus annotatedVirus : reportableViruses) {
            oncogenic.add("virus");

            ActionabilityKey keyVirus = ImmutableActionabilityKey.builder()
                    .match(annotatedVirus.interpretation() != null ? annotatedVirus.interpretation() : Strings.EMPTY)
                    .type(TypeAlteration.POSITIVE)
                    .build();
            ActionabilityEntry entry = actionabilityMap.get(keyVirus);
            if (entry != null && entry.condition() == Condition.ALWAYS) {
                conclusion.put(conclusion.size(), "- " + annotatedVirus.interpretation() + " " + entry.conclusion());
                actionable.add("virus");
            }
        }
    }

    public static void generateHrdConclusion(@NotNull Map<Integer, String> conclusion, @NotNull ChordData chordAnalysis,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap, @NotNull Set<String> oncogenic,
            @NotNull Set<String> actionable, @NotNull Set<String> HRD) {
        if (chordAnalysis.hrStatus() == ChordStatus.HR_DEFICIENT) {
            ActionabilityKey keyHRD = ImmutableActionabilityKey.builder().match("HRD").type(TypeAlteration.POSITIVE).build();
            ActionabilityEntry entry = actionabilityMap.get(keyHRD);
            if (entry != null && entry.condition() == Condition.ALWAYS) {
                if (HRD.size() == 0) {
                    ActionabilityKey keyNoHRD =
                            ImmutableActionabilityKey.builder().match("NO_HRD_CAUSE").type(TypeAlteration.NO_HRD_CAUSE).build();
                    ActionabilityEntry entryNoHRd = actionabilityMap.get(keyNoHRD);
                    if (entryNoHRd != null && entry.condition() == Condition.OTHER) {
                        conclusion.put(conclusion.size(),
                                "- " + "HRD (" + chordAnalysis.hrdValue() + ") " + entry.conclusion() + entryNoHRd.conclusion());
                    }
                }
                conclusion.put(conclusion.size(),
                        "- " + "HRD (" + DOUBLE_DECIMAL_FORMAT.format(chordAnalysis.hrdValue()) + ") " + entry.conclusion());

                actionable.add("HRD");
                oncogenic.add("HRD");
            }
        }
    }

    public static void generateMSIConclusion(@NotNull Map<Integer, String> conclusion, @NotNull MicrosatelliteStatus microsatelliteStatus,
            double microsatelliteMb, @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap, @NotNull Set<String> oncogenic,
            @NotNull Set<String> actionable) {
        if (microsatelliteStatus == MicrosatelliteStatus.MSI) {
            ActionabilityKey keyMSI = ImmutableActionabilityKey.builder().match("MSI").type(TypeAlteration.POSITIVE).build();
            ActionabilityEntry entry = actionabilityMap.get(keyMSI);
            if (entry != null && entry.condition() == Condition.ALWAYS) {
                conclusion.put(conclusion.size(),
                        "- " + "MSI (" + DOUBLE_DECIMAL_FORMAT.format(microsatelliteMb) + ") " + entry.conclusion());
                actionable.add("MSI");
                oncogenic.add("MSI");
            }
        }
    }

    public static void generateTMLConclusion(@NotNull Map<Integer, String> conclusion, @NotNull TumorMutationalStatus tumorMutationalStatus,
            int tumorMutationalLoad, @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap, @NotNull Set<String> oncogenic,
            @NotNull Set<String> actionable) {
        if (tumorMutationalStatus == TumorMutationalStatus.HIGH) {
            ActionabilityKey keyTML = ImmutableActionabilityKey.builder().match("High-TML").type(TypeAlteration.POSITIVE).build();
            ActionabilityEntry entry = actionabilityMap.get(keyTML);
            if (entry != null && entry.condition() == Condition.ALWAYS) {
                conclusion.put(conclusion.size(), "- " + "TML (" + tumorMutationalLoad + ") " + entry.conclusion());
                actionable.add("TML");
                oncogenic.add("TML");
            }
        }
    }

    public static void generateTMBConclusion(@NotNull Map<Integer, String> conclusion, double tumorMutationalBurden,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap, @NotNull Set<String> oncogenic,
            @NotNull Set<String> actionable) {
        if (tumorMutationalBurden >= TMB_CUTOFF) {
            ActionabilityKey keyTMB = ImmutableActionabilityKey.builder().match("High-TMB").type(TypeAlteration.POSITIVE).build();
            ActionabilityEntry entry = actionabilityMap.get(keyTMB);
            if (entry != null && entry.condition() == Condition.ALWAYS) {
                conclusion.put(conclusion.size(), "- " + "TMB (" + tumorMutationalBurden + ") " + entry.conclusion());
                actionable.add("TMB");
                oncogenic.add("TMB");
            }
        }
    }

    public static void genertatePurityConclusion(@NotNull Map<Integer, String> conclusion, double purity, boolean hasRelaiblePurity,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap) {
        if (!hasRelaiblePurity) {
            ActionabilityKey keyReliable =
                    ImmutableActionabilityKey.builder().match("PURITY_UNRELIABLE").type(TypeAlteration.PURITY_UNRELIABLE).build();

            ActionabilityEntry entryReliable = actionabilityMap.get(keyReliable);
            if (entryReliable != null && entryReliable.condition() == Condition.OTHER) {
                conclusion.put(conclusion.size(), "- " + entryReliable.conclusion());
            }
        } else if (purity < PURITY_CUTOFF) {
            ActionabilityKey keyPurity = ImmutableActionabilityKey.builder().match("PURITY").type(TypeAlteration.PURITY).build();

            ActionabilityEntry entry = actionabilityMap.get(keyPurity);
            if (entry != null && entry.condition() == Condition.OTHER) {
                conclusion.put(conclusion.size(), "- " + entry.conclusion().replace("XX%", purity + "%"));
            }
        }
    }

    public static void generateTotalResults(@NotNull Map<Integer, String> conclusion,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap, @NotNull Set<String> oncogenic,
            @NotNull Set<String> actionable) {
        if (oncogenic.size() == 0) {
            ActionabilityKey keyOncogenic =
                    ImmutableActionabilityKey.builder().match("NO_ONCOGENIC").type(TypeAlteration.NO_ONCOGENIC).build();

            ActionabilityEntry entry = actionabilityMap.get(keyOncogenic);
            if (entry != null && entry.condition() == Condition.OTHER) {
                conclusion.put(conclusion.size(), "- " + entry.conclusion());
            }
        } else if (actionable.size() == 0) {
            ActionabilityKey keyActionable =
                    ImmutableActionabilityKey.builder().match("NO_ACTIONABLE").type(TypeAlteration.NO_ACTIONABLE).build();
            ActionabilityEntry entry = actionabilityMap.get(keyActionable);
            if (entry != null && entry.condition() == Condition.OTHER) {
                conclusion.put(conclusion.size(), "- " + entry.conclusion());
            }
        }
    }

    public static void generateFindings(@NotNull Map<Integer, String> conclusion,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap) {
        ActionabilityKey keyOncogenic = ImmutableActionabilityKey.builder().match("FINDINGS").type(TypeAlteration.FINDINGS).build();

        ActionabilityEntry entry = actionabilityMap.get(keyOncogenic);
        if (entry != null && entry.condition() == Condition.OTHER) {
            conclusion.put(conclusion.size(), "- " + entry.conclusion());
        }
    }

    @NotNull
    public static DecimalFormat decimalFormat(@NotNull String format) {
        // To make sure every decimal format uses a dot as separator rather than a comma.
        return new DecimalFormat(format, DecimalFormatSymbols.getInstance(Locale.ENGLISH));
    }
}