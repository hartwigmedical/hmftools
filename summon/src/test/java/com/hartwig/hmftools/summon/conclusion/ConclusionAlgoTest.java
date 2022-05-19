package com.hartwig.hmftools.summon.conclusion;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.Map;
import java.util.Set;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.summon.actionability.ActionabilityEntry;
import com.hartwig.hmftools.summon.actionability.ActionabilityKey;
import com.hartwig.hmftools.summon.actionability.ImmutableActionabilityEntry;
import com.hartwig.hmftools.summon.actionability.ImmutableActionabilityKey;
import com.hartwig.hmftools.summon.actionability.Type;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ConclusionAlgoTest {

    @Test
    public void canGenerateConclusion(){
    }

    @Test
    public void canGenerateConclusionString(){
    }

    @Test
    public void canGenerateActionabilityMap(){
    }

    @Test
    public void canGenerateDriverGenesMap(){
    }

    @Test
    public void canGenerateCUPPAConclusion(){

    }

    @Test
    public void canGenerateSomaticConclusion(){

    }

    @Test
    public void canGenerateGermlineConclusion(){

    }

    @Test
    public void canGenerateCNVConclusion(){

    }

    @Test
    public void canGenerateFusionConclusion(){

    }

    @Test
    public void canGenerateHomozygousDisruptionConclusion(){

    }

    @Test
    public void canGenerateVirusConclusion(){

    }

    @Test
    public void canGenerateHrdConclusion(){

    }

    @Test
    public void canGenerateMSIConclusion(){

    }

    @Test
    public void canGenerateTMLConclusion(){

    }

    @Test
    public void canGenerateTMBConclusion(){

    }

    @Test
    public void canGenertatePurityConclusionBelow(){
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        ActionabilityKey key = ImmutableActionabilityKey.builder().gene("purity").type(Type.PURITY).build();
        ActionabilityEntry entry = ImmutableActionabilityEntry.builder().gene("purity").type(Type.PURITY).onlyHighDriver(false).conclusion("low purity (XX%)").build();
        actionabilityMap.put(key, entry);
        ConclusionAlgo.genertatePurityConclusion(conclusion, 0.1, actionabilityMap);
        assertEquals(conclusion.get(0), "- low purity (0.1%)");
    }

    @Test
    public void canGenertatePurityConclusionAbove(){
        Map<Integer, String> conclusion = Maps.newHashMap();
        Map<ActionabilityKey, ActionabilityEntry> actionabilityMap = Maps.newHashMap();
        ActionabilityKey key = ImmutableActionabilityKey.builder().gene("purity").type(Type.PURITY).build();
        ActionabilityEntry entry = ImmutableActionabilityEntry.builder().gene("purity").type(Type.PURITY).onlyHighDriver(false).conclusion("low purity (XX%)").build();
        actionabilityMap.put(key, entry);
        ConclusionAlgo.genertatePurityConclusion(conclusion, 0.3, actionabilityMap);
        assertNull(conclusion.get(0));
    }

    @Test
    public void canGenerateTotalResults(){

    }

    @Test
    public void canGenerateFindings(){

    }
}