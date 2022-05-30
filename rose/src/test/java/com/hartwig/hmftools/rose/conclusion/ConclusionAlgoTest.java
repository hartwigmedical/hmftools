package com.hartwig.hmftools.rose.conclusion;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.chord.ChordTestFactory;
import com.hartwig.hmftools.common.chord.ImmutableChordAnalysis;
import com.hartwig.hmftools.common.cuppa.ImmutableMolecularTissueOrigin;
import com.hartwig.hmftools.common.cuppa.MolecularTissueOrigin;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogTestFactory;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.linx.ImmutableReportableHomozygousDisruption;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.linx.ReportableHomozygousDisruption;
import com.hartwig.hmftools.common.purple.PurpleTestFactory;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ImmutableReportableGainLoss;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.sv.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.test.SomaticVariantTestFactory;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
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
    public void canGenerateConclusionString() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        conclusion.put(0, "fusion");
        conclusion.put(1, "variant");
        conclusion.put(2, "variant");
        conclusion.put(3, "variant");
        conclusion.put(4, "loss");
        conclusion.put(5, "amplification");
        assertEquals(ConclusionAlgo.generateConslusionString(conclusion),
                "fusion <enter> variant <enter> variant <enter> variant <enter> loss <enter> amplification <enter> ");
    }

    @Test
    public void canGenerateCUPPAConclusion() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "CUPPA", TypeAlteration.CUPPA, "CUPPA", Condition.OTHER, "CUPPA");

        MolecularTissueOrigin molecularTissueOrigin =
                ImmutableMolecularTissueOrigin.builder().plotPath(Strings.EMPTY).conclusion("Melanoma").build();
        ConclusionAlgo.generateCUPPAConclusion(conclusion, molecularTissueOrigin, actionabilityMap);
        assertEquals(conclusion.get(0), "- CUPPA");
    }

    @Test
    public void canGenerateCUPPAConclusionInconclusive() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();

        actionabilityMap = testActionabilityMap(actionabilityMap,
                "CUPPA_INCONCLUSIVE",
                TypeAlteration.CUPPA_INCONCLUSIVE,
                "CUPPA_INCONCLUSIVE",
                Condition.OTHER,
                "results inconclusive");

        MolecularTissueOrigin molecularTissueOrigin =
                ImmutableMolecularTissueOrigin.builder().plotPath(Strings.EMPTY).conclusion("results inconclusive").build();
        ConclusionAlgo.generateCUPPAConclusion(conclusion, molecularTissueOrigin, actionabilityMap);
        assertEquals(conclusion.get(0), "- results inconclusive");
    }

    @Test
    public void canGenerateVariantsConclusion() {
        List<ReportableVariant> reportableVariants = canGenerateVariants();
        Map<String, DriverGene> driverGenesMap = Maps.newHashMap();
        driverGenesMap.put("CHEK2", createDriverGene("CHEK2", DriverCategory.TSG));
        driverGenesMap.put("APC", createDriverGene("APC", ONCO));
        driverGenesMap.put("BRCA2", createDriverGene("BRCA2", DriverCategory.TSG));
        driverGenesMap.put("BRCA1", createDriverGene("BRCA1", DriverCategory.TSG));

        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "CHEK2", TypeAlteration.INACTIVATION, "CHEK2", Condition.ONLY_HIGH, "CHEK2");
        actionabilityMap = testActionabilityMap(actionabilityMap, "APC", TypeAlteration.ACTIVATING_MUTATION, "APC", Condition.ALWAYS_NO_ACTIONABLE, "APC");
        actionabilityMap = testActionabilityMap(actionabilityMap, "BRCA2", TypeAlteration.INACTIVATION, "BRCA2", Condition.ONLY_HIGH, "BRCA2");
        actionabilityMap = testActionabilityMap(actionabilityMap, "BRCA1", TypeAlteration.INACTIVATION, "BRCA1", Condition.ONLY_HIGH, "BRCA1");
        actionabilityMap = testActionabilityMap(actionabilityMap, "germline", TypeAlteration.GERMLINE, "germline", Condition.ONLY_HIGH, "germline");
        actionabilityMap = testActionabilityMap(actionabilityMap, "NOT_BIALLELIC", TypeAlteration.NOT_BIALLELIC, "NOT_BIALLELIC", Condition.OTHER, "not biallelic");

        ConclusionAlgo.generateVariantConclusion(conclusion,
                reportableVariants,
                actionabilityMap,
                driverGenesMap,
                Sets.newHashSet(),
                Sets.newHashSet(),
                Sets.newHashSet());
        assertEquals(conclusion.size(), 3);
        assertEquals(conclusion.get(0), "- CHEK2(p.?) CHEK2 not biallelic");
        assertEquals(conclusion.get(1), "- APC(p.?) APC");
        assertEquals(conclusion.get(2), "- BRCA2(p.?) BRCA2");
        assertNull(conclusion.get(3));
    }

    @Test
    public void canGenerateCNVConclusion() {
        List<ReportableGainLoss> gainLosse = gainloss();
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();

        actionabilityMap = testActionabilityMap(actionabilityMap, "BRAF", TypeAlteration.AMPLIFICATION, "BRAF", Condition.ALWAYS, "BRAF");
        actionabilityMap = testActionabilityMap(actionabilityMap, "KRAS", TypeAlteration.AMPLIFICATION, "KRAS", Condition.ALWAYS, "KRAS");
        actionabilityMap = testActionabilityMap(actionabilityMap, "CDKN2A", TypeAlteration.LOSS, "CDKN2A", Condition.ALWAYS, "CDKN2A");
        actionabilityMap = testActionabilityMap(actionabilityMap, "EGFR", TypeAlteration.LOSS, "EGFR", Condition.ALWAYS, "EGFR");

        ConclusionAlgo.generateCNVConclusion(conclusion, gainLosse, actionabilityMap, Sets.newHashSet(), Sets.newHashSet());
        assertEquals(conclusion.size(), 4);
        assertEquals(conclusion.get(0), "- BRAF BRAF");
        assertEquals(conclusion.get(1), "- KRAS KRAS");
        assertEquals(conclusion.get(2), "- CDKN2A CDKN2A");
        assertEquals(conclusion.get(3), "- EGFR EGFR");
    }

    @Test
    public void canGenerateFusionConclusion() {
        List<LinxFusion> fusions = createFusion();
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "BRAF", TypeAlteration.INTERNAL_DELETION, "BRAF", Condition.ALWAYS, "BRAF");
        actionabilityMap = testActionabilityMap(actionabilityMap, "MET", TypeAlteration.FUSION, "MET", Condition.ALWAYS, "MET");
        actionabilityMap = testActionabilityMap(actionabilityMap, "EGFR", TypeAlteration.KINASE_DOMAIN_DUPLICATION, "EGFR", Condition.ALWAYS, "EGFR");

        ConclusionAlgo.generateFusionConclusion(conclusion, fusions, actionabilityMap, Sets.newHashSet(), Sets.newHashSet());
        assertEquals(conclusion.size(), 4);
        assertEquals(conclusion.get(0), "- BRAF-BRAF BRAF");
        assertEquals(conclusion.get(1), "- CAV2-MET MET");
        assertEquals(conclusion.get(2), "- EGFR-EGFR EGFR");
        assertEquals(conclusion.get(3), "- EGFR-EGFR EGFR");

    }

    @Test
    public void canGenerateHomozygousDisruptionConclusion() {
        List<ReportableHomozygousDisruption> homozygousDisruptions = Lists.newArrayList(createHomozygousDisruption("PTEN"));
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "PTEN", TypeAlteration.INACTIVATION, "PTEN", Condition.ALWAYS, "PTEN");

        ConclusionAlgo.generateHomozygousDisruptionConclusion(conclusion,
                homozygousDisruptions,
                actionabilityMap,
                Sets.newHashSet(),
                Sets.newHashSet());

        assertEquals(conclusion.size(), 1);
        assertEquals(conclusion.get(0), "- PTEN PTEN");
    }

    @Test
    public void canGenerateVirusConclusion() {
        List<AnnotatedVirus> annotatedVirus = annotatedVirus();
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();

        actionabilityMap = testActionabilityMap(actionabilityMap, "EBV", TypeAlteration.POSITIVE, "EBV", Condition.ALWAYS, "EBV");
        actionabilityMap = testActionabilityMap(actionabilityMap, "HPV", TypeAlteration.POSITIVE, "HPV", Condition.ALWAYS, "HPV");

        ConclusionAlgo.generateVirusConclusion(conclusion, annotatedVirus, actionabilityMap, Sets.newHashSet(), Sets.newHashSet());
        assertEquals(conclusion.size(), 2);
        assertEquals(conclusion.get(0), "- EBV EBV");
        assertEquals(conclusion.get(1), "- HPV HPV");
    }

    @Test
    public void canGenerateHrdConclusion() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "HRD", TypeAlteration.POSITIVE, "HRD", Condition.ALWAYS, "HRD");

        ChordAnalysis analysis = ImmutableChordAnalysis.builder()
                .from(ChordTestFactory.createMinimalTestChordAnalysis())
                .hrdValue(0.8)
                .hrStatus(ChordStatus.HR_DEFICIENT)
                .build();
        Set<String> HRD = Sets.newHashSet();
        HRD.add("BRCA1");
        ConclusionAlgo.generateHrdConclusion(conclusion, analysis, actionabilityMap, Sets.newHashSet(), Sets.newHashSet(), HRD);
        assertEquals(conclusion.get(0), "- HRD(0.8) HRD");
    }

    @Test
    public void canGenerateHrpConclusion() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "HRD", TypeAlteration.POSITIVE, "HRD", Condition.ALWAYS, "HRD");

        ChordAnalysis analysis = ImmutableChordAnalysis.builder()
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
        assertNull(conclusion.get(0));
    }

    @Test
    public void canGenerateMSIConclusion() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "MSI", TypeAlteration.POSITIVE, "MSI", Condition.ALWAYS, "MSI");

        ConclusionAlgo.generateMSIConclusion(conclusion,
                MicrosatelliteStatus.MSI,
                4.5,
                actionabilityMap,
                Sets.newHashSet(),
                Sets.newHashSet());
        assertEquals(conclusion.get(0), "- MSI(4.5)MSI");
    }

    @Test
    public void canGenerateMSSConclusion() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "MSI", TypeAlteration.POSITIVE, "MSI", Condition.ALWAYS, "MSI");

        ConclusionAlgo.generateMSIConclusion(conclusion,
                MicrosatelliteStatus.MSS,
                3.2,
                actionabilityMap,
                Sets.newHashSet(),
                Sets.newHashSet());
        assertNull(conclusion.get(0));
    }

    @Test
    public void canGenerateTMLHighConclusion() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "High-TML", TypeAlteration.POSITIVE, "High-TML", Condition.ALWAYS, "TML");

        ConclusionAlgo.generateTMLConclusion(conclusion,
                TumorMutationalStatus.HIGH,
                200,
                actionabilityMap,
                Sets.newHashSet(),
                Sets.newHashSet());
        assertEquals(conclusion.get(0), "- TML(200) TML");
    }

    @Test
    public void canGenerateTMLLowConclusion() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "High-TML", TypeAlteration.POSITIVE, "High-TML", Condition.ALWAYS, "TML");

        ConclusionAlgo.generateTMLConclusion(conclusion,
                TumorMutationalStatus.LOW,
                100,
                actionabilityMap,
                Sets.newHashSet(),
                Sets.newHashSet());
        assertNull(conclusion.get(0));
    }

    @Test
    public void canGenerateTMBHighConclusion() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "High-TMB", TypeAlteration.POSITIVE, "High-TMB", Condition.ALWAYS, "TMB");

        ConclusionAlgo.generateTMBConclusion(conclusion, 15, actionabilityMap, Sets.newHashSet(), Sets.newHashSet());
        assertEquals(conclusion.get(0), "- TMB( 15.0)TMB");
    }

    @Test
    public void canGenerateTMBLowConclusion() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap, "High-TMB", TypeAlteration.POSITIVE, "High-TMB", Condition.ALWAYS, "TMB");

        ConclusionAlgo.generateTMBConclusion(conclusion, 9, actionabilityMap, Sets.newHashSet(), Sets.newHashSet());
        assertNull(conclusion.get(0));
    }

    @Test
    public void canGenertatePurityConclusionBelow() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap =
                testActionabilityMap(actionabilityMap, "PURITY", TypeAlteration.PURITY, "PURITY", Condition.OTHER, "low purity (XX%)");

        ConclusionAlgo.genertatePurityConclusion(conclusion, 0.1, true, actionabilityMap);
        assertEquals(conclusion.get(0), "- low purity (0.1%)");
    }

    @Test
    public void canGenertatePurityConclusionAbove() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap =
                testActionabilityMap(actionabilityMap, "PURITY", TypeAlteration.PURITY, "PURITY", Condition.OTHER, "low purity (XX%)");

        ConclusionAlgo.genertatePurityConclusion(conclusion, 0.3, true, actionabilityMap);
        assertNull(conclusion.get(0));
    }

    @Test
    public void canGenertatePurityConclusionReliable() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap,
                "PURITY_UNRELIABLE",
                TypeAlteration.PURITY_UNRELIABLE,
                "PURITY_UNRELIABLE",
                Condition.OTHER,
                "unreliable");

        ConclusionAlgo.genertatePurityConclusion(conclusion, 0.3, false, actionabilityMap);
        assertEquals(conclusion.get(0), "- unreliable");
    }

    @Test
    public void canGenerateTotalResultsOncogenic() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap,
                "NO_ONCOGENIC",
                TypeAlteration.NO_ONCOGENIC,
                "NO_ONCOGENIC",
                Condition.OTHER,
                "no_oncogenic");

        ConclusionAlgo.generateTotalResults(conclusion, actionabilityMap, Sets.newHashSet(), Sets.newHashSet());
        assertEquals(conclusion.get(0), "- no_oncogenic");
    }

    @Test
    public void canGenerateTotalResultsActionable() {
        Set<String> oncogenic = Sets.newHashSet();
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap,
                "NO_ACTIONABLE",
                TypeAlteration.NO_ACTIONABLE,
                "NO_ACTIONABLE",
                Condition.OTHER,
                "no_actionable");

        oncogenic.add("fusion");
        ConclusionAlgo.generateTotalResults(conclusion, actionabilityMap, oncogenic, Sets.newHashSet());
        assertEquals(conclusion.get(0), "- no_actionable");

    }

    @Test
    public void canGenerateFindings() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        actionabilityMap = testActionabilityMap(actionabilityMap,
                "FINDINGS",
                TypeAlteration.FINDINGS,
                "FINDINGS",
                Condition.OTHER,
                "findings");

        ConclusionAlgo.generateFindings(conclusion, actionabilityMap);
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
                .canonicalHgvsProteinImpact("p.?")
                .biallelic(true)
                .build();

        SomaticVariant variant2 = SomaticVariantTestFactory.builder()
                .reported(true)
                .gene("BRCA2")
                .canonicalTranscript("transcript1")
                .canonicalHgvsProteinImpact("p.?")
                .biallelic(true)
                .build();

        SomaticVariant variant3 = SomaticVariantTestFactory.builder()
                .reported(true)
                .gene("BRCA1")
                .canonicalTranscript("transcript1")
                .canonicalHgvsProteinImpact("p.?")
                .biallelic(true)
                .build();

        List<ReportableVariant> reportableSomatic =
                ReportableVariantFactory.toReportableSomaticVariants(Lists.newArrayList(variant1, variant2, variant3),
                        Lists.newArrayList(DriverCatalogTestFactory.createCanonicalSomaticMutationEntryForGene("APC",
                                        0.4,
                                        "transcript1",
                                        ONCO),
                                DriverCatalogTestFactory.createCanonicalSomaticMutationEntryForGene("BRCA2", 0.9, "transcript1", TSG),
                                DriverCatalogTestFactory.createCanonicalSomaticMutationEntryForGene("BRCA1", 0.7, "transcript1", TSG)));

        SomaticVariant variant4 = SomaticVariantTestFactory.builder()
                .reported(true)
                .gene("CHEK2")
                .canonicalTranscript("transcript1")
                .canonicalHgvsProteinImpact("p.?")
                .biallelic(false)
                .build();

        List<ReportableVariant> reportableGermline = ReportableVariantFactory.toReportableGermlineVariants(Lists.newArrayList(variant4),
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
                .build();
    }

    @NotNull
    public List<ReportableGainLoss> gainloss() {
        ReportableGainLoss gainLoss1 = ImmutableReportableGainLoss.builder()
                .from(PurpleTestFactory.createReportableGainLoss("BRAF", CopyNumberInterpretation.FULL_GAIN))
                .build();
        ReportableGainLoss gainLoss2 = ImmutableReportableGainLoss.builder()
                .from(PurpleTestFactory.createReportableGainLoss("KRAS", CopyNumberInterpretation.PARTIAL_GAIN))
                .build();
        ReportableGainLoss gainLoss3 = ImmutableReportableGainLoss.builder()
                .from(PurpleTestFactory.createReportableGainLoss("CDKN2A", CopyNumberInterpretation.FULL_LOSS))
                .build();
        ReportableGainLoss gainLoss4 = ImmutableReportableGainLoss.builder()
                .from(PurpleTestFactory.createReportableGainLoss("EGFR", CopyNumberInterpretation.PARTIAL_LOSS))
                .build();
        return Lists.newArrayList(gainLoss1, gainLoss2, gainLoss3, gainLoss4);
    }

    @NotNull
    public List<AnnotatedVirus> annotatedVirus() {
        List<AnnotatedVirus> annotatedVirus = Lists.newArrayList();
        AnnotatedVirus virus1 = VirusTestFactory.testAnnotatedVirusBuilder().interpretation("EBV").build();
        AnnotatedVirus virus2 = VirusTestFactory.testAnnotatedVirusBuilder().interpretation("HPV").build();
        annotatedVirus.add(virus1);
        annotatedVirus.add(virus2);
        return annotatedVirus;
    }

    @NotNull
    private static ReportableHomozygousDisruption createHomozygousDisruption(@NotNull String gene) {
        return ImmutableReportableHomozygousDisruption.builder()
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