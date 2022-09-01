package com.hartwig.hmftools.rose.conclusion;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

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
import com.hartwig.hmftools.common.rose.ActionabilityConclusion;
import com.hartwig.hmftools.common.rose.ImmutableActionabilityConclusion;
import com.hartwig.hmftools.common.variant.DriverInterpretation;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantFactory;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.purple.TumorMutationalStatus;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;
import com.hartwig.hmftools.rose.RoseData;
import com.hartwig.hmftools.rose.actionability.ActionabilityEntry;
import com.hartwig.hmftools.rose.actionability.ActionabilityKey;
import com.hartwig.hmftools.rose.actionability.Condition;
import com.hartwig.hmftools.rose.actionability.ImmutableActionabilityKey;
import com.hartwig.hmftools.rose.actionability.TypeAlteration;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ConclusionAlgo {
    private static final Logger LOGGER = LogManager.getLogger(ConclusionAlgo.class);

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
        List<String> conclusion = Lists.newArrayList();
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

        generateCUPPAConclusion(conclusion, roseData.molecularTissueOrigin(), actionabilityMap);
        generateVariantConclusion(conclusion,
                reportableVariants,
                actionabilityMap,
                driverGenesMap,
                oncogenic,
                actionable,
                HRD,
                roseData.chord());
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
        generatePurityConclusion(conclusion, roseData.purple().purity(), roseData.purple().hasReliablePurity(), actionabilityMap);
        generateFindings(conclusion, actionabilityMap);

        return ImmutableActionabilityConclusion.builder().conclusion(conclusion).build();
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

    public static void generateCUPPAConclusion(@NotNull List<String> conclusion, MolecularTissueOrigin molecularTissueOrigin,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap) {

        if (molecularTissueOrigin.conclusion().contains("results inconclusive")) {
            ActionabilityKey keyCuppaInconclusice =
                    ImmutableActionabilityKey.builder().match("CUPPA_INCONCLUSIVE").type(TypeAlteration.CUPPA_INCONCLUSIVE).build();

            ActionabilityEntry entry = actionabilityMap.get(keyCuppaInconclusice);
            if (entry != null && entry.condition() == Condition.OTHER) {
                conclusion.add("- " + entry.conclusion());
            }
        } else {
            ActionabilityKey keyCuppa = ImmutableActionabilityKey.builder().match("CUPPA").type(TypeAlteration.CUPPA).build();

            ActionabilityEntry entry = actionabilityMap.get(keyCuppa);
            if (entry != null && entry.condition() == Condition.OTHER) {
                conclusion.add("- " + entry.conclusion() + " " + molecularTissueOrigin.conclusion());
            }
        }
    }

    public static void generateVariantConclusion(@NotNull List<String> conclusion, @NotNull List<ReportableVariant> reportableVariants,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap, @NotNull Map<String, DriverGene> driverGenesMap,
            @NotNull Set<String> oncogenic, @NotNull Set<String> actionable, @NotNull Set<String> HRD, @NotNull ChordData chordAnalysis) {

        Map<String, List<VariantKey>> variantKeyList = Maps.newHashMap();

        for (ReportableVariant reportableVariant : reportableVariants) {
            List<VariantKey> variantKeys = Lists.newArrayList();
            VariantKey variantKey = ImmutableVariantKey.builder()
                    .gene(reportableVariant.gene())
                    .variantAnnotation(EventGenerator.variantEvent(reportableVariant))
                    .driverInterpretation(reportableVariant.driverLikelihoodInterpretation())
                    .bialleic(reportableVariant.biallelic())
                    .build();

            if (variantKeyList.containsKey(reportableVariant.gene())) {
                variantKeys.addAll(variantKeyList.get(reportableVariant.gene()));
                variantKeys.add(variantKey);
                variantKeyList.put(reportableVariant.gene(), variantKeys);
            } else {
                variantKeys.add(variantKey);
                variantKeyList.put(reportableVariant.gene(), variantKeys);
            }
        }

        for (Map.Entry<String, List<VariantKey>> keyMap : variantKeyList.entrySet()) {
            boolean HRDgene = false;
            TypeAlteration alteration = TypeAlteration.UNKNOWN;

            Set<Boolean> biallelic = Sets.newHashSet();
            StringJoiner variantMerging = new StringJoiner(",");
            for (VariantKey key : keyMap.getValue()) {
                if (HRD_GENES.contains(keyMap.getKey())) {
                    HRD.add(keyMap.getKey());
                    HRDgene = true;
                }
                oncogenic.add("variant");

                if (keyMap.getKey().equals("KRAS") && key.variantAnnotation().equals("p.Gly12Cys")) {
                    alteration = TypeAlteration.ACTIVATING_MUTATION_KRAS_G12C;
                } else if (driverGenesMap.get(key.gene()).likelihoodType().equals(DriverCategory.ONCO)) {
                    alteration = TypeAlteration.ACTIVATING_MUTATION;
                } else if (driverGenesMap.get(key.gene()).likelihoodType().equals(DriverCategory.TSG)) {
                    alteration = TypeAlteration.INACTIVATION;
                }

                variantMerging.add(key.variantAnnotation());
                biallelic.add(key.bialleic());
            }

            ActionabilityKey keySomaticVariant = ImmutableActionabilityKey.builder().match(keyMap.getKey()).type(alteration).build();
            ActionabilityEntry entry = actionabilityMap.get(keySomaticVariant);
            if (entry != null) {
                if ((keyMap.getValue().get(0).driverInterpretation() == DriverInterpretation.HIGH
                        && entry.condition() == Condition.ONLY_HIGH) || entry.condition() == Condition.ALWAYS_NO_ACTIONABLE) {
                    if (entry.condition() == Condition.ONLY_HIGH) {
                        actionable.add("variant");
                    }
                    boolean biallelicBoolean = false;
                    if (biallelic.size() == 1) {
                        biallelicBoolean = biallelic.stream().findFirst().get();
                    }

                    if (driverGenesMap.get(keyMap.getKey()).likelihoodType().equals(DriverCategory.TSG) && !biallelicBoolean) {
                        ActionabilityKey keyBiallelic =
                                ImmutableActionabilityKey.builder().match("NOT_BIALLELIC").type(TypeAlteration.NOT_BIALLELIC).build();
                        ActionabilityEntry entryBiallelic = actionabilityMap.get(keyBiallelic);
                        if (entryBiallelic.condition() == Condition.OTHER) {
                            conclusion.add("- " + keyMap.getKey() + " (" + variantMerging + ") " + entry.conclusion() + " "
                                    + entryBiallelic.conclusion());
                        }
                    } else {
                        conclusion.add("- " + keyMap.getKey() + " (" + variantMerging + ") " + entry.conclusion());
                    }
                } else if (HRDgene && chordAnalysis.hrStatus() == ChordStatus.HR_DEFICIENT) {
                    conclusion.add("- " + keyMap.getKey() + " (" + variantMerging + ") " + entry.conclusion());
                    actionable.add("variant");
                }
            }
        }
    }

    public static void generateCNVConclusion(@NotNull List<String> conclusion, @NotNull List<GainLoss> reportableGainLosses,
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
                    conclusion.add("- " + gainLoss.gene() + " " + entry.conclusion());
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
                    conclusion.add("- " + gainLoss.gene() + " " + entry.conclusion());
                    actionable.add("CNV");
                }
            }
        }
    }

    public static void generateFusionConclusion(@NotNull List<String> conclusion, @NotNull List<LinxFusion> reportableFusions,
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
                    conclusion.add("- " + fusion.geneStart() + " - " + fusion.geneEnd() + " (" + fusion.geneContextStart() + " - "
                            + fusion.geneContextEnd() + ") " + entry.conclusion());
                    actionable.add("fusion");
                }
            } else if (fusion.reportedType().equals(KnownFusionType.EXON_DEL_DUP.toString())) {
                ActionabilityKey keyFusion =
                        ImmutableActionabilityKey.builder().match(fusion.geneStart()).type(TypeAlteration.INTERNAL_DELETION).build();
                ActionabilityEntry entry = actionabilityMap.get(keyFusion);
                if (entry != null && entry.condition() == Condition.ALWAYS) {
                    conclusion.add("- " + fusion.geneStart() + " - " + fusion.geneEnd() + " (" + fusion.geneContextStart() + " - "
                            + fusion.geneContextEnd() + ") " + entry.conclusion());
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
                    conclusion.add("- " + fusion.geneStart() + " - " + fusion.geneEnd() + " (" + fusion.geneContextStart() + " - "
                            + fusion.geneContextEnd() + ") " + entryStart.conclusion());
                    actionable.add("fusion");
                } else if (entryEnd != null && entryEnd.condition() == Condition.ALWAYS) {
                    conclusion.add("- " + fusion.geneStart() + " - " + fusion.geneEnd() + " (" + fusion.geneContextStart() + " - "
                            + fusion.geneContextEnd() + ") " + entryEnd.conclusion());
                    actionable.add("fusion");
                }
            }
        }
    }

    public static void generateHomozygousDisruptionConclusion(@NotNull List<String> conclusion,
            @NotNull List<HomozygousDisruption> homozygousDisruptions, @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap,
            @NotNull Set<String> oncogenic, @NotNull Set<String> actionable) {

        for (HomozygousDisruption homozygousDisruption : homozygousDisruptions) {
            oncogenic.add("homozygousDisruption");

            ActionabilityKey keyHomozygousDisruption =
                    ImmutableActionabilityKey.builder().match(homozygousDisruption.gene()).type(TypeAlteration.INACTIVATION).build();
            ActionabilityEntry entry = actionabilityMap.get(keyHomozygousDisruption);
            if (entry != null && entry.condition() == Condition.ALWAYS) {
                conclusion.add("- " + homozygousDisruption.gene() + " " + entry.conclusion());
                actionable.add("homozygousDisruption");
            }
        }
    }

    public static void generateVirusConclusion(@NotNull List<String> conclusion, @NotNull List<AnnotatedVirus> reportableViruses,
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
                conclusion.add("- " + annotatedVirus.interpretation() + " " + entry.conclusion());
                actionable.add("virus");
            } else if (entry == null) {
                if (annotatedVirus.interpretation() != null && (annotatedVirus.virusDriverLikelihoodType() == VirusLikelihoodType.LOW
                        || annotatedVirus.virusDriverLikelihoodType() == VirusLikelihoodType.HIGH)) {
                    conclusion.add("- " + annotatedVirus.interpretation() + " positive");
                }
            }
        }
    }

    public static void generateHrdConclusion(@NotNull List<String> conclusion, @NotNull ChordData chordAnalysis,
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
                        conclusion.add("- " + "HRD (" + chordAnalysis.hrdValue() + ") " + entry.conclusion() + entryNoHRd.conclusion());
                    }
                }
                conclusion.add("- " + "HRD (" + DOUBLE_DECIMAL_FORMAT.format(chordAnalysis.hrdValue()) + ") " + entry.conclusion());

                actionable.add("HRD");
                oncogenic.add("HRD");
            }
        }
    }

    public static void generateMSIConclusion(@NotNull List<String> conclusion, @NotNull MicrosatelliteStatus microsatelliteStatus,
            double microsatelliteMb, @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap, @NotNull Set<String> oncogenic,
            @NotNull Set<String> actionable) {

        if (microsatelliteStatus == MicrosatelliteStatus.MSI) {
            ActionabilityKey keyMSI = ImmutableActionabilityKey.builder().match("MSI").type(TypeAlteration.POSITIVE).build();
            ActionabilityEntry entry = actionabilityMap.get(keyMSI);
            if (entry != null && entry.condition() == Condition.ALWAYS) {
                conclusion.add("- " + "MSI (" + DOUBLE_DECIMAL_FORMAT.format(microsatelliteMb) + ") " + entry.conclusion());
                actionable.add("MSI");
                oncogenic.add("MSI");
            }
        }
    }

    public static void generateTMLConclusion(@NotNull List<String> conclusion, @NotNull TumorMutationalStatus tumorMutationalStatus,
            int tumorMutationalLoad, @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap, @NotNull Set<String> oncogenic,
            @NotNull Set<String> actionable) {

        if (tumorMutationalStatus == TumorMutationalStatus.HIGH) {
            ActionabilityKey keyTML = ImmutableActionabilityKey.builder().match("High-TML").type(TypeAlteration.POSITIVE).build();
            ActionabilityEntry entry = actionabilityMap.get(keyTML);
            if (entry != null && entry.condition() == Condition.ALWAYS) {
                conclusion.add("- " + "TML (" + tumorMutationalLoad + ") " + entry.conclusion());
                actionable.add("TML");
                oncogenic.add("TML");
            }
        }
    }

    public static void generateTMBConclusion(@NotNull List<String> conclusion, double tumorMutationalBurden,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap, @NotNull Set<String> oncogenic,
            @NotNull Set<String> actionable) {

        if (tumorMutationalBurden >= TMB_CUTOFF) {
            ActionabilityKey keyTMB = ImmutableActionabilityKey.builder().match("High-TMB").type(TypeAlteration.POSITIVE).build();
            ActionabilityEntry entry = actionabilityMap.get(keyTMB);
            if (entry != null && entry.condition() == Condition.ALWAYS) {
                conclusion.add("- " + "TMB (" + tumorMutationalBurden + ") " + entry.conclusion());
                actionable.add("TMB");
                oncogenic.add("TMB");
            }
        }
    }

    public static void generatePurityConclusion(@NotNull List<String> conclusion, double purity, boolean hasRelaiblePurity,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap) {
        if (!hasRelaiblePurity) {
            ActionabilityKey keyReliable =
                    ImmutableActionabilityKey.builder().match("PURITY_UNRELIABLE").type(TypeAlteration.PURITY_UNRELIABLE).build();

            ActionabilityEntry entryReliable = actionabilityMap.get(keyReliable);
            if (entryReliable != null && entryReliable.condition() == Condition.OTHER) {
                conclusion.add("- " + entryReliable.conclusion());
            }
        } else if (purity < PURITY_CUTOFF) {
            ActionabilityKey keyPurity = ImmutableActionabilityKey.builder().match("PURITY").type(TypeAlteration.PURITY).build();

            ActionabilityEntry entry = actionabilityMap.get(keyPurity);
            if (entry != null && entry.condition() == Condition.OTHER) {
                conclusion.add("- " + entry.conclusion().replace("XX%", purity + "%"));
            }
        }
    }

    public static void generateTotalResults(@NotNull List<String> conclusion,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap, @NotNull Set<String> oncogenic,
            @NotNull Set<String> actionable) {

        if (oncogenic.size() == 0) {
            ActionabilityKey keyOncogenic =
                    ImmutableActionabilityKey.builder().match("NO_ONCOGENIC").type(TypeAlteration.NO_ONCOGENIC).build();

            ActionabilityEntry entry = actionabilityMap.get(keyOncogenic);
            if (entry != null && entry.condition() == Condition.OTHER) {
                conclusion.add("- " + entry.conclusion());
            }
        } else if (actionable.size() == 0) {
            ActionabilityKey keyActionable =
                    ImmutableActionabilityKey.builder().match("NO_ACTIONABLE").type(TypeAlteration.NO_ACTIONABLE).build();
            ActionabilityEntry entry = actionabilityMap.get(keyActionable);
            if (entry != null && entry.condition() == Condition.OTHER) {
                conclusion.add("- " + entry.conclusion());
            }
        }
    }

    public static void generateFindings(@NotNull List<String> conclusion,
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap) {
        ActionabilityKey keyOncogenic = ImmutableActionabilityKey.builder().match("FINDINGS").type(TypeAlteration.FINDINGS).build();

        ActionabilityEntry entry = actionabilityMap.get(keyOncogenic);
        if (entry != null && entry.condition() == Condition.OTHER) {
            conclusion.add("- " + entry.conclusion());
        }
    }

    @NotNull
    public static DecimalFormat decimalFormat(@NotNull String format) {
        // To make sure every decimal format uses a dot as separator rather than a comma.
        return new DecimalFormat(format, DecimalFormatSymbols.getInstance(Locale.ENGLISH));
    }
}