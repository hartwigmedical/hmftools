package com.hartwig.hmftools.rose.conclusion;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.chord.ChordTestFactory;
import com.hartwig.hmftools.common.chord.ImmutableChordData;
import com.hartwig.hmftools.common.cuppa.interpretation.CuppaPrediction;
import com.hartwig.hmftools.common.cuppa.interpretation.ImmutableCuppaPrediction;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogTestFactory;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.linx.HomozygousDisruption;
import com.hartwig.hmftools.common.linx.ImmutableHomozygousDisruption;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.purple.loader.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.loader.GainLoss;
import com.hartwig.hmftools.common.purple.loader.GainLossTestFactory;
import com.hartwig.hmftools.common.purple.loader.ImmutableGainLoss;
import com.hartwig.hmftools.common.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.test.SomaticVariantTestFactory;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.purple.TumorMutationalStatus;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;
import com.hartwig.hmftools.common.virus.VirusTestFactory;
import com.hartwig.hmftools.rose.actionability.ActionabilityEntry;
import com.hartwig.hmftools.rose.actionability.ActionabilityKey;
import com.hartwig.hmftools.rose.actionability.Condition;
import com.hartwig.hmftools.rose.actionability.ImmutableActionabilityEntry;
import com.hartwig.hmftools.rose.actionability.ImmutableActionabilityKey;
import com.hartwig.hmftools.rose.actionability.TypeAlteration;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ConclusionAlgoTest {

    @Test
    public void canGenerateCUPPAConclusion() {
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap,
                "CUPPA",
                TypeAlteration.CUPPA,
                "CUPPA",
                Condition.OTHER,
                "Molecular Tissue of Origin classifier: XXXX.");

        CuppaPrediction cuppaPrediction = ImmutableCuppaPrediction.builder().cancerType("Melanoma").likelihood(0.996).build();

        ConclusionAlgo.generateCUPPAConclusion(conclusion, cuppaPrediction, actionabilityMap);

        assertEquals(1, conclusion.size());
        assertEquals(conclusion.get(0), "- Molecular Tissue of Origin classifier: Melanoma (likelihood: 99.6%).");
    }

    @Test
    public void canGenerateCUPPAConclusionInconclusive() {
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();

        actionabilityMap = testActionabilityMap(actionabilityMap,
                "CUPPA_INCONCLUSIVE",
                TypeAlteration.CUPPA_INCONCLUSIVE,
                "CUPPA_INCONCLUSIVE",
                Condition.OTHER,
                "Molecular Tissue of Origin classifier: Inconclusive (highest likelihood: xxx - xx%).");

        CuppaPrediction cuppaPrediction = ImmutableCuppaPrediction.builder().cancerType("Melanoma").likelihood(0.45).build();

        ConclusionAlgo.generateCUPPAConclusion(conclusion, cuppaPrediction, actionabilityMap);

        assertEquals(1, conclusion.size());
        assertEquals(conclusion.get(0), "- Molecular Tissue of Origin classifier: Inconclusive.");
    }

    @Test
    public void canGenerateCUPPAConclusionInconclusiveWithLocation() {
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();

        actionabilityMap = testActionabilityMap(actionabilityMap,
                "CUPPA_INCONCLUSIVE",
                TypeAlteration.CUPPA_INCONCLUSIVE,
                "CUPPA_INCONCLUSIVE",
                Condition.OTHER,
                "Molecular Tissue of Origin classifier: Inconclusive (highest likelihood: xxx - xx%).");

        CuppaPrediction cuppaPrediction = ImmutableCuppaPrediction.builder().cancerType("Melanoma").likelihood(0.601).build();

        ConclusionAlgo.generateCUPPAConclusion(conclusion, cuppaPrediction, actionabilityMap);

        assertEquals(1, conclusion.size());
        assertEquals(conclusion.get(0), "- Molecular Tissue of Origin classifier: Inconclusive (highest likelihood: Melanoma-60.1%).");
    }

    @Test
    public void canGenerateVariantsConclusion() {
        List<ReportableVariant> reportableVariants = canGenerateVariants();
        Map<String, DriverGene> driverGenesMap = Maps.newHashMap();
        driverGenesMap.put("CHEK2", createDriverGene("CHEK2", DriverCategory.TSG));
        driverGenesMap.put("APC", createDriverGene("APC", ONCO));
        driverGenesMap.put("BRCA2", createDriverGene("BRCA2", DriverCategory.TSG));
        driverGenesMap.put("BRCA1", createDriverGene("BRCA1", DriverCategory.TSG));

        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap =
                testActionabilityMap(actionabilityMap, "CHEK2", TypeAlteration.INACTIVATION, "CHEK2", Condition.ONLY_HIGH, "CHEK2");
        actionabilityMap = testActionabilityMap(actionabilityMap,
                "APC",
                TypeAlteration.ACTIVATING_MUTATION,
                "APC",
                Condition.ALWAYS_NO_ACTIONABLE,
                "APC");
        actionabilityMap =
                testActionabilityMap(actionabilityMap, "BRCA2", TypeAlteration.INACTIVATION, "BRCA2", Condition.ONLY_HIGH, "BRCA2");
        actionabilityMap =
                testActionabilityMap(actionabilityMap, "BRCA1", TypeAlteration.INACTIVATION, "BRCA1", Condition.ONLY_HIGH, "BRCA1");
        actionabilityMap =
                testActionabilityMap(actionabilityMap, "germline", TypeAlteration.GERMLINE, "germline", Condition.ONLY_HIGH, "germline");
        actionabilityMap = testActionabilityMap(actionabilityMap,
                "NOT_BIALLELIC",
                TypeAlteration.NOT_BIALLELIC,
                "NOT_BIALLELIC",
                Condition.OTHER,
                "not biallelic");

        ChordData analysis = ImmutableChordData.builder()
                .from(ChordTestFactory.createMinimalTestChordAnalysis())
                .hrdValue(0.8)
                .hrStatus(ChordStatus.HR_DEFICIENT)
                .build();

        ConclusionAlgo.generateVariantConclusion(conclusion,
                reportableVariants,
                actionabilityMap,
                driverGenesMap,
                Sets.newHashSet(),
                Sets.newHashSet(),
                Sets.newHashSet(),
                analysis);
        assertEquals(4, conclusion.size());
        assertEquals(conclusion.get(0), "- APC (p.Val600Arg) APC");
        assertEquals(conclusion.get(1), "- CHEK2 (c.123A>C splice) CHEK2 not biallelic");
        assertEquals(conclusion.get(2), "- BRCA1 (p.Val600Arg,p.Val602Arg) BRCA1");
        assertEquals(conclusion.get(3), "- BRCA2 (c.1235A>C splice) BRCA2");
    }

    @Test
    public void canGenerateCNVConclusion() {
        List<GainLoss> gainLosse = gainloss();
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();

        actionabilityMap = testActionabilityMap(actionabilityMap, "BRAF", TypeAlteration.AMPLIFICATION, "BRAF", Condition.ALWAYS, "BRAF");
        actionabilityMap = testActionabilityMap(actionabilityMap, "KRAS", TypeAlteration.AMPLIFICATION, "KRAS", Condition.ALWAYS, "KRAS");
        actionabilityMap = testActionabilityMap(actionabilityMap, "CDKN2A", TypeAlteration.LOSS, "CDKN2A", Condition.ALWAYS, "CDKN2A");
        actionabilityMap = testActionabilityMap(actionabilityMap, "EGFR", TypeAlteration.LOSS, "EGFR", Condition.ALWAYS, "EGFR");

        ConclusionAlgo.generateCNVConclusion(conclusion, gainLosse, actionabilityMap, Sets.newHashSet(), Sets.newHashSet());
        assertEquals(4, conclusion.size());
        assertEquals(conclusion.get(0), "- BRAF (copies: 4) BRAF");
        assertEquals(conclusion.get(1), "- KRAS (copies: 8) KRAS");
        assertEquals(conclusion.get(2), "- CDKN2A (copies: 0) CDKN2A");
        assertEquals(conclusion.get(3), "- EGFR (copies: 0) EGFR");
    }

    @Test
    public void canGenerateFusionConclusion() {
        List<LinxFusion> fusions = createFusion();
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap =
                testActionabilityMap(actionabilityMap, "BRAF", TypeAlteration.INTERNAL_DELETION, "BRAF", Condition.ALWAYS, "BRAF");
        actionabilityMap = testActionabilityMap(actionabilityMap, "MET", TypeAlteration.FUSION, "MET", Condition.ALWAYS, "MET");
        actionabilityMap =
                testActionabilityMap(actionabilityMap, "EGFR", TypeAlteration.KINASE_DOMAIN_DUPLICATION, "EGFR", Condition.ALWAYS, "EGFR");

        ConclusionAlgo.generateFusionConclusion(conclusion, fusions, actionabilityMap, Sets.newHashSet(), Sets.newHashSet());
        assertEquals(4, conclusion.size());
        assertEquals(conclusion.get(0), "- BRAF - BRAF ( - ) BRAF");
        assertEquals(conclusion.get(1), "- CAV2 - MET ( - ) MET");
        assertEquals(conclusion.get(2), "- EGFR - EGFR ( - ) EGFR");
        assertEquals(conclusion.get(3), "- EGFR - EGFR ( - ) EGFR");
    }

    @Test
    public void canGenerateHomozygousDisruptionConclusion() {
        List<HomozygousDisruption> homozygousDisruptions = Lists.newArrayList(createHomozygousDisruption("PTEN"));
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "PTEN", TypeAlteration.INACTIVATION, "PTEN", Condition.ALWAYS, "PTEN");

        ConclusionAlgo.generateHomozygousDisruptionConclusion(conclusion,
                homozygousDisruptions,
                actionabilityMap,
                Sets.newHashSet(),
                Sets.newHashSet());

        assertEquals(1, conclusion.size());
        assertEquals(conclusion.get(0), "- PTEN PTEN");
    }

    @Test
    public void canGenerateVirusConclusion() {
        List<AnnotatedVirus> annotatedVirus = annotatedVirus();
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();

        actionabilityMap = testActionabilityMap(actionabilityMap, "EBV", TypeAlteration.POSITIVE, "EBV", Condition.ALWAYS, "EBV");
        actionabilityMap = testActionabilityMap(actionabilityMap, "HPV", TypeAlteration.POSITIVE, "HPV", Condition.ALWAYS, "HPV");

        ConclusionAlgo.generateVirusConclusion(conclusion, annotatedVirus, actionabilityMap, Sets.newHashSet(), Sets.newHashSet());
        assertEquals(3, conclusion.size());
        assertEquals(conclusion.get(0), "- EBV EBV");
        assertEquals(conclusion.get(1), "- HPV HPV");
        assertEquals(conclusion.get(2), "- MCV positive");
    }

    @Test
    public void canGenerateHrdConclusion() {
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "HRD", TypeAlteration.POSITIVE, "HRD", Condition.ALWAYS, "HRD");

        ChordData analysis = ImmutableChordData.builder()
                .from(ChordTestFactory.createMinimalTestChordAnalysis())
                .hrdValue(0.8)
                .hrStatus(ChordStatus.HR_DEFICIENT)
                .build();
        Set<String> HRD = Sets.newHashSet();
        HRD.add("BRCA1");
        ConclusionAlgo.generateHrdConclusion(conclusion, analysis, actionabilityMap, Sets.newHashSet(), Sets.newHashSet(), HRD);

        assertEquals(1, conclusion.size());
        assertEquals(conclusion.get(0), "- HRD (0.8) HRD");
    }

    @Test
    public void canGenerateHrpConclusion() {
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "HRD", TypeAlteration.POSITIVE, "HRD", Condition.ALWAYS, "HRD");

        ChordData analysis = ImmutableChordData.builder()
                .from(ChordTestFactory.createMinimalTestChordAnalysis())
                .hrdValue(0.4)
                .hrStatus(ChordStatus.HR_PROFICIENT)
                .build();
        ConclusionAlgo.generateHrdConclusion(conclusion,
                analysis,
                actionabilityMap,
                Sets.newHashSet(),
                Sets.newHashSet(),
                Sets.newHashSet());

        assertEquals(0, conclusion.size());
    }

    @Test
    public void canGenerateMSIConclusion() {
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "MSI", TypeAlteration.POSITIVE, "MSI", Condition.ALWAYS, "MSI");

        ConclusionAlgo.generateMSIConclusion(conclusion,
                MicrosatelliteStatus.MSI,
                4.5,
                actionabilityMap,
                Sets.newHashSet(),
                Sets.newHashSet());

        assertEquals(1, conclusion.size());
        assertEquals(conclusion.get(0), "- MSI (4.5) MSI");
    }

    @Test
    public void canGenerateMSSConclusion() {
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "MSI", TypeAlteration.POSITIVE, "MSI", Condition.ALWAYS, "MSI");

        ConclusionAlgo.generateMSIConclusion(conclusion,
                MicrosatelliteStatus.MSS,
                3.2,
                actionabilityMap,
                Sets.newHashSet(),
                Sets.newHashSet());
        assertEquals(0, conclusion.size());
    }

    @Test
    public void canGenerateTMLHighConclusion() {
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "High-TML", TypeAlteration.POSITIVE, "High-TML", Condition.ALWAYS, "TML");

        ConclusionAlgo.generateTMLConclusion(conclusion,
                TumorMutationalStatus.HIGH,
                200,
                actionabilityMap,
                Sets.newHashSet(),
                Sets.newHashSet());

        assertEquals(1, conclusion.size());
        assertEquals(conclusion.get(0), "- TML (200) TML");
    }

    @Test
    public void canGenerateTMLLowConclusion() {
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "High-TML", TypeAlteration.POSITIVE, "High-TML", Condition.ALWAYS, "TML");

        ConclusionAlgo.generateTMLConclusion(conclusion,
                TumorMutationalStatus.LOW,
                100,
                actionabilityMap,
                Sets.newHashSet(),
                Sets.newHashSet());

        assertEquals(0, conclusion.size());
    }

    @Test
    public void canGenerateTMBHighConclusion() {
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "High-TMB", TypeAlteration.POSITIVE, "High-TMB", Condition.ALWAYS, "TMB");

        ConclusionAlgo.generateTMBConclusion(conclusion, 15, actionabilityMap, Sets.newHashSet(), Sets.newHashSet());

        assertEquals(1, conclusion.size());
        assertEquals(conclusion.get(0), "- TMB (15.0) TMB");
    }

    @Test
    public void canGenerateTMBLowConclusion() {
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "High-TMB", TypeAlteration.POSITIVE, "High-TMB", Condition.ALWAYS, "TMB");

        ConclusionAlgo.generateTMBConclusion(conclusion, 9, actionabilityMap, Sets.newHashSet(), Sets.newHashSet());
        assertEquals(0, conclusion.size());
    }

    @Test
    public void canGenertatePurityConclusionBelow() {
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap =
                testActionabilityMap(actionabilityMap, "PURITY", TypeAlteration.PURITY, "PURITY", Condition.OTHER, "low purity (XX%)");

        ConclusionAlgo.generatePurityConclusion(conclusion, 0.16, true, actionabilityMap);

        assertEquals(1, conclusion.size());
        assertEquals(conclusion.get(0), "- low purity (16%)");
    }

    @Test
    public void canGenertatePurityConclusionAbove() {
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap =
                testActionabilityMap(actionabilityMap, "PURITY", TypeAlteration.PURITY, "PURITY", Condition.OTHER, "low purity (XX%)");

        ConclusionAlgo.generatePurityConclusion(conclusion, 0.3, true, actionabilityMap);
        assertEquals(0, conclusion.size());
    }

    @Test
    public void canGenertatePurityConclusionReliable() {
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap,
                "PURITY_UNRELIABLE",
                TypeAlteration.PURITY_UNRELIABLE,
                "PURITY_UNRELIABLE",
                Condition.OTHER,
                "unreliable");

        ConclusionAlgo.generatePurityConclusion(conclusion, 0.3, false, actionabilityMap);

        assertEquals(1, conclusion.size());
        assertEquals(conclusion.get(0), "- unreliable");
    }

    @Test
    public void canGenerateTotalResultsOncogenic() {
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap,
                "NO_ONCOGENIC",
                TypeAlteration.NO_ONCOGENIC,
                "NO_ONCOGENIC",
                Condition.OTHER,
                "no_oncogenic");

        ConclusionAlgo.generateTotalResults(conclusion, actionabilityMap, Sets.newHashSet(), Sets.newHashSet());

        assertEquals(1, conclusion.size());
        assertEquals(conclusion.get(0), "- no_oncogenic");
    }

    @Test
    public void canGenerateTotalResultsActionable() {
        Set<String> oncogenic = Sets.newHashSet();
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap,
                "NO_ACTIONABLE",
                TypeAlteration.NO_ACTIONABLE,
                "NO_ACTIONABLE",
                Condition.OTHER,
                "no_actionable");

        oncogenic.add("fusion");
        ConclusionAlgo.generateTotalResults(conclusion, actionabilityMap, oncogenic, Sets.newHashSet());

        assertEquals(1, conclusion.size());
        assertEquals(conclusion.get(0), "- no_actionable");

    }

    @Test
    public void canGenerateFindings() {
        List<String> conclusion = Lists.newArrayList();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap =
                testActionabilityMap(actionabilityMap, "FINDINGS", TypeAlteration.FINDINGS, "FINDINGS", Condition.OTHER, "findings");

        ConclusionAlgo.generateFindings(conclusion, actionabilityMap);

        assertEquals(1, conclusion.size());
        assertEquals(conclusion.get(0), "- findings");
    }

    @NotNull
    public Map<ActionabilityKey, ActionabilityEntry> testActionabilityMap(
            @NotNull Map<ActionabilityKey, ActionabilityEntry> actionabilityMap, @NotNull String gene,
            @NotNull TypeAlteration typeAlteration, @NotNull String match, @NotNull Condition condition, @NotNull String conclusion) {
        ActionabilityKey key = ImmutableActionabilityKey.builder().match(gene).type(typeAlteration).build();
        ActionabilityEntry entry =
                ImmutableActionabilityEntry.builder().match(match).type(typeAlteration).condition(condition).conclusion(conclusion).build();
        actionabilityMap.put(key, entry);
        return actionabilityMap;
    }

    @NotNull
    public List<ReportableVariant> canGenerateVariants() {
        SomaticVariant variant1 = SomaticVariantTestFactory.builder()
                .reported(true)
                .gene("APC")
                .canonicalTranscript("transcript1")
                .canonicalHgvsProteinImpact("p.Val600Arg")
                .canonicalHgvsCodingImpact("c.123A>C")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .biallelic(true)
                .build();

        SomaticVariant variant2 = SomaticVariantTestFactory.builder()
                .reported(true)
                .gene("BRCA2")
                .canonicalTranscript("transcript1")
                .canonicalHgvsProteinImpact("p.?")
                .canonicalHgvsCodingImpact("c.1235A>C")
                .canonicalCodingEffect(CodingEffect.SPLICE)
                .biallelic(true)
                .build();

        SomaticVariant variant3 = SomaticVariantTestFactory.builder()
                .reported(true)
                .gene("BRCA1")
                .canonicalTranscript("transcript1")
                .canonicalHgvsProteinImpact("p.Val600Arg")
                .canonicalHgvsCodingImpact("c.123A>C")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .biallelic(true)
                .build();

        SomaticVariant variant4 = SomaticVariantTestFactory.builder()
                .reported(true)
                .gene("BRCA1")
                .canonicalTranscript("transcript1")
                .canonicalHgvsProteinImpact("p.Val602Arg")
                .canonicalHgvsCodingImpact("c.124A>C")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .biallelic(true)
                .build();

        List<ReportableVariant> reportableSomatic =
                ReportableVariantFactory.toReportableSomaticVariants(Lists.newArrayList(variant1, variant2, variant3, variant4),
                        Lists.newArrayList(DriverCatalogTestFactory.createCanonicalSomaticMutationEntryForGene("APC",
                                        0.4,
                                        "transcript1",
                                        ONCO),
                                DriverCatalogTestFactory.createCanonicalSomaticMutationEntryForGene("BRCA2", 0.9, "transcript1", TSG),
                                DriverCatalogTestFactory.createCanonicalSomaticMutationEntryForGene("BRCA1", 0.82, "transcript1", TSG)));

        SomaticVariant variant5 = SomaticVariantTestFactory.builder()
                .reported(true)
                .gene("CHEK2")
                .canonicalTranscript("transcript1")
                .canonicalHgvsProteinImpact("")
                .canonicalHgvsCodingImpact("c.123A>C")
                .canonicalCodingEffect(CodingEffect.SPLICE)
                .biallelic(false)
                .build();

        List<ReportableVariant> reportableGermline = ReportableVariantFactory.toReportableGermlineVariants(Lists.newArrayList(variant5),
                Lists.newArrayList(DriverCatalogTestFactory.createCanonicalGermlineMutationEntryForGene("CHEK2",
                        0.85,
                        "transcript1",
                        TSG)));

        return ReportableVariantFactory.mergeVariantLists(reportableGermline, reportableSomatic);
    }

    public static DriverGene createDriverGene(@NotNull String name, @NotNull DriverCategory likelihoodMethod) {
        return ImmutableDriverGene.builder()
                .gene(name)
                .reportMissenseAndInframe(false)
                .reportNonsenseAndFrameshift(false)
                .reportSplice(false)
                .reportDeletion(false)
                .reportDisruption(true)
                .reportAmplification(false)
                .reportSomaticHotspot(false)
                .reportGermlineVariant(DriverGeneGermlineReporting.NONE)
                .reportGermlineHotspot(DriverGeneGermlineReporting.NONE)
                .likelihoodType(likelihoodMethod)
                .reportGermlineDisruption(true)
                .reportPGX(false)
                .build();
    }

    @NotNull
    public List<GainLoss> gainloss() {
        GainLoss gainLoss1 = ImmutableGainLoss.builder()
                .from(GainLossTestFactory.createGainLoss("BRAF", CopyNumberInterpretation.FULL_GAIN))
                .minCopies(4)
                .maxCopies(4)
                .build();
        GainLoss gainLoss2 = ImmutableGainLoss.builder()
                .from(GainLossTestFactory.createGainLoss("KRAS", CopyNumberInterpretation.PARTIAL_GAIN))
                .minCopies(3)
                .maxCopies(8)
                .build();
        GainLoss gainLoss3 = ImmutableGainLoss.builder()
                .from(GainLossTestFactory.createGainLoss("CDKN2A", CopyNumberInterpretation.FULL_LOSS))
                .minCopies(0)
                .maxCopies(0)
                .build();
        GainLoss gainLoss4 = ImmutableGainLoss.builder()
                .from(GainLossTestFactory.createGainLoss("EGFR", CopyNumberInterpretation.PARTIAL_LOSS))
                .minCopies(0)
                .maxCopies(3)
                .build();
        return Lists.newArrayList(gainLoss1, gainLoss2, gainLoss3, gainLoss4);
    }

    @NotNull
    public List<AnnotatedVirus> annotatedVirus() {
        List<AnnotatedVirus> annotatedVirus = Lists.newArrayList();
        AnnotatedVirus virus1 = VirusTestFactory.testAnnotatedVirusBuilder().interpretation("EBV").build();
        AnnotatedVirus virus2 = VirusTestFactory.testAnnotatedVirusBuilder().interpretation("HPV").build();
        AnnotatedVirus virus3 = VirusTestFactory.testAnnotatedVirusBuilder()
                .interpretation("MCV")
                .virusDriverLikelihoodType(VirusLikelihoodType.HIGH)
                .build();
        annotatedVirus.add(virus1);
        annotatedVirus.add(virus2);
        annotatedVirus.add(virus3);
        return annotatedVirus;
    }

    @NotNull
    private static HomozygousDisruption createHomozygousDisruption(@NotNull String gene) {
        return ImmutableHomozygousDisruption.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(gene)
                .transcript("123")
                .isCanonical(true)
                .build();
    }

    @NotNull
    private static List<LinxFusion> createFusion() {
        List<LinxFusion> fusion = Lists.newArrayList();
        fusion.add(linxFusionBuilder("BRAF", "BRAF", true).reportedType(KnownFusionType.EXON_DEL_DUP.toString()).name("BRAF-BRAF").build());
        fusion.add(linxFusionBuilder("CAV2", "MET", true).reportedType(KnownFusionType.KNOWN_PAIR.toString()).name("CAV2-MET").build());
        fusion.add(linxFusionBuilder("EGFR", "EGFR", true).fusedExonUp(25)
                .fusedExonDown(14)
                .reportedType(KnownFusionType.EXON_DEL_DUP.toString())
                .name("EGFR-EGFR")
                .build());
        fusion.add(linxFusionBuilder("EGFR", "EGFR", true).fusedExonUp(26)
                .fusedExonDown(18)
                .reportedType(KnownFusionType.EXON_DEL_DUP.toString())
                .name("EGFR-EGFR")
                .build());
        fusion.add(linxFusionBuilder("EGFR", "EGFR", true).fusedExonUp(15)
                .fusedExonDown(23)
                .reportedType(KnownFusionType.EXON_DEL_DUP.toString())
                .name("EGFR-EGFR")
                .build());
        return fusion;
    }

    @NotNull
    private static ImmutableLinxFusion.Builder linxFusionBuilder(@NotNull String geneStart, @NotNull String geneEnd, boolean reported) {
        return ImmutableLinxFusion.builder()
                .from(LinxTestFactory.createMinimalTestFusion())
                .geneStart(geneStart)
                .geneEnd(geneEnd)
                .reported(reported);
    }
}