package com.hartwig.hmftools.summon.conclusion;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.Map;
import java.util.Set;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.chord.ChordTestFactory;
import com.hartwig.hmftools.common.chord.ImmutableChordAnalysis;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.summon.actionability.ActionabilityEntry;
import com.hartwig.hmftools.summon.actionability.ActionabilityKey;
import com.hartwig.hmftools.summon.actionability.ImmutableActionabilityEntry;
import com.hartwig.hmftools.summon.actionability.ImmutableActionabilityKey;
import com.hartwig.hmftools.summon.actionability.Type;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ConclusionAlgoTest {

    @Test
    public void canGenerateConclusion() {
    }

    @Test
    public void canGenerateConclusionString() {
    }

    @Test
    public void canGenerateActionabilityMap() {
    }

    @Test
    public void canGenerateDriverGenesMap() {
    }

    @Test
    public void canGenerateCUPPAConclusion() {

    }

    @Test
    public void canGenerateSomaticConclusion() {

    }

    @Test
    public void canGenerateGermlineConclusion() {

    }

    @Test
    public void canGenerateCNVConclusion() {

    }

    @Test
    public void canGenerateFusionConclusion() {

    }

    @Test
    public void canGenerateHomozygousDisruptionConclusion() {

    }

    @Test
    public void canGenerateVirusConclusion() {

    }

    @Test
    public void canGenerateHrdConclusion() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        ActionabilityKey key = ImmutableActionabilityKey.builder().gene("HRD").type(Type.POSITIVE).build();
        ActionabilityEntry entry = ImmutableActionabilityEntry.builder()
                .gene("HRD")
                .type(Type.POSITIVE)
                .onlyHighDriver(false)
                .conclusion("HRD")
                .build();
        actionabilityMap.put(key, entry);

        ChordAnalysis analysis = ImmutableChordAnalysis.builder()
                .from(ChordTestFactory.createMinimalTestChordAnalysis())
                .hrdValue(0.8)
                .hrStatus(ChordStatus.HR_DEFICIENT)
                .build();
        ConclusionAlgo.generateHrdConclusion(conclusion, analysis, actionabilityMap);
        assertEquals(conclusion.get(0), "- HRD(0.8) HRD");
    }

    @Test
    public void canGenerateHrpConclusion() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        ActionabilityKey key = ImmutableActionabilityKey.builder().gene("HRD").type(Type.POSITIVE).build();
        ActionabilityEntry entry = ImmutableActionabilityEntry.builder()
                .gene("HRD")
                .type(Type.POSITIVE)
                .onlyHighDriver(false)
                .conclusion("HRD")
                .build();
        actionabilityMap.put(key, entry);

        ChordAnalysis analysis = ImmutableChordAnalysis.builder()
                .from(ChordTestFactory.createMinimalTestChordAnalysis())
                .hrdValue(0.4)
                .hrStatus(ChordStatus.HR_PROFICIENT)
                .build();
        ConclusionAlgo.generateHrdConclusion(conclusion, analysis, actionabilityMap);
        assertNull(conclusion.get(0));
    }

    @Test
    public void canGenerateMSIConclusion() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        ActionabilityKey key = ImmutableActionabilityKey.builder().gene("MSI").type(Type.POSITIVE).build();
        ActionabilityEntry entry = ImmutableActionabilityEntry.builder()
                .gene("MSI")
                .type(Type.POSITIVE)
                .onlyHighDriver(false)
                .conclusion("MSI")
                .build();
        actionabilityMap.put(key, entry);
        ConclusionAlgo.generateMSIConclusion(conclusion, MicrosatelliteStatus.MSI, 4.5, actionabilityMap);
        assertEquals(conclusion.get(0), "- MSI(4.5)MSI");
    }

    @Test
    public void canGenerateMSSConclusion() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        ActionabilityKey key = ImmutableActionabilityKey.builder().gene("MSI").type(Type.POSITIVE).build();
        ActionabilityEntry entry = ImmutableActionabilityEntry.builder()
                .gene("MSI")
                .type(Type.POSITIVE)
                .onlyHighDriver(false)
                .conclusion("MSI")
                .build();
        actionabilityMap.put(key, entry);
        ConclusionAlgo.generateMSIConclusion(conclusion, MicrosatelliteStatus.MSS, 3.2, actionabilityMap);
        assertNull(conclusion.get(0));
    }

    @Test
    public void canGenerateTMLHighConclusion() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        ActionabilityKey key = ImmutableActionabilityKey.builder().gene("High-TML").type(Type.POSITIVE).build();
        ActionabilityEntry entry = ImmutableActionabilityEntry.builder()
                .gene("High-TML")
                .type(Type.POSITIVE)
                .onlyHighDriver(false)
                .conclusion("TML")
                .build();
        actionabilityMap.put(key, entry);
        ConclusionAlgo.generateTMLConclusion(conclusion, TumorMutationalStatus.HIGH, 200, actionabilityMap);
        assertEquals(conclusion.get(0), "- TML(200) TML");
    }

    @Test
    public void canGenerateTMLLowConclusion() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        ActionabilityKey key = ImmutableActionabilityKey.builder().gene("High-TML").type(Type.POSITIVE).build();
        ActionabilityEntry entry = ImmutableActionabilityEntry.builder()
                .gene("High-TML")
                .type(Type.POSITIVE)
                .onlyHighDriver(false)
                .conclusion("TML")
                .build();
        actionabilityMap.put(key, entry);
        ConclusionAlgo.generateTMLConclusion(conclusion, TumorMutationalStatus.LOW, 100, actionabilityMap);
        assertNull(conclusion.get(0));
    }

    @Test
    public void canGenerateTMBHighConclusion() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        ActionabilityKey key = ImmutableActionabilityKey.builder().gene("High-TMB").type(Type.POSITIVE).build();
        ActionabilityEntry entry = ImmutableActionabilityEntry.builder()
                .gene("High-TMB")
                .type(Type.POSITIVE)
                .onlyHighDriver(false)
                .conclusion("TMB")
                .build();
        actionabilityMap.put(key, entry);
        ConclusionAlgo.generateTMBConclusion(conclusion, 15, actionabilityMap);
        assertEquals(conclusion.get(0), "- TMB( 15.0)TMB");
    }

    @Test
    public void canGenerateTMBLowConclusion() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        ActionabilityKey key = ImmutableActionabilityKey.builder().gene("High-TMB").type(Type.POSITIVE).build();
        ActionabilityEntry entry = ImmutableActionabilityEntry.builder()
                .gene("High-TMB")
                .type(Type.POSITIVE)
                .onlyHighDriver(false)
                .conclusion("TMB")
                .build();
        actionabilityMap.put(key, entry);
        ConclusionAlgo.generateTMBConclusion(conclusion, 15, actionabilityMap);
        assertNull(conclusion.get(0));
    }

    @Test
    public void canGenertatePurityConclusionBelow() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        ActionabilityKey key = ImmutableActionabilityKey.builder().gene("purity").type(Type.PURITY).build();
        ActionabilityEntry entry = ImmutableActionabilityEntry.builder()
                .gene("purity")
                .type(Type.PURITY)
                .onlyHighDriver(false)
                .conclusion("low purity (XX%)")
                .build();
        actionabilityMap.put(key, entry);
        ConclusionAlgo.genertatePurityConclusion(conclusion, 0.1, actionabilityMap);
        assertEquals(conclusion.get(0), "- low purity (0.1%)");
    }

    @Test
    public void canGenertatePurityConclusionAbove() {
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        ActionabilityKey key = ImmutableActionabilityKey.builder().gene("purity").type(Type.PURITY).build();
        ActionabilityEntry entry = ImmutableActionabilityEntry.builder()
                .gene("purity")
                .type(Type.PURITY)
                .onlyHighDriver(false)
                .conclusion("low purity (XX%)")
                .build();
        actionabilityMap.put(key, entry);
        ConclusionAlgo.genertatePurityConclusion(conclusion, 0.3, actionabilityMap);
        assertNull(conclusion.get(0));
    }

    @Test
    public void canGenerateTotalResults() {

    }

    @Test
    public void canGenerateFindings() {

    }
}