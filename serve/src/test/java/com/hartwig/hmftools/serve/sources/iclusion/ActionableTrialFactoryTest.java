package com.hartwig.hmftools.serve.sources.iclusion;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTumorLocation;
import com.hartwig.hmftools.iclusion.datamodel.ImmutableIclusionTumorLocation;
import com.hartwig.hmftools.serve.curation.DoidLookupTestFactory;

import org.junit.Test;

public class ActionableTrialFactoryTest {

    @Test
    public void canCreateActionableTrials() {
        String location1 = "loc1";
        String loc1Doid1 = "loc1Doid1";
        String loc1Doid2 = "loc1Doid2";
        String location2 = "loc2";
        String loc2Doid1 = "loc2Doid2";
        String treatment = "trial";

        IclusionTumorLocation loc1 =
                ImmutableIclusionTumorLocation.builder().primaryTumorLocation(location1).addDoids(loc1Doid1).addDoids(loc1Doid2).build();
        IclusionTumorLocation loc2 = ImmutableIclusionTumorLocation.builder().primaryTumorLocation(location2).addDoids(loc2Doid1).build();
        IclusionTrial trial = IclusionTestFactory.trialWithTumors(treatment, Lists.newArrayList(loc1, loc2));

        ActionableTrialFactory factory = new ActionableTrialFactory(DoidLookupTestFactory.dummy());
        List<ActionableTrial> actionableTrials = factory.toActionableTrials(trial);
        assertEquals(3, actionableTrials.size());
        assertEquals(treatment, actionableTrials.get(0).treatment());
        assertEquals(location1, actionableTrials.get(0).cancerType());
        assertEquals(loc1Doid1, actionableTrials.get(0).doid());

        assertEquals(treatment, actionableTrials.get(1).treatment());
        assertEquals(location1, actionableTrials.get(1).cancerType());
        assertEquals(loc1Doid2, actionableTrials.get(1).doid());

        assertEquals(treatment, actionableTrials.get(2).treatment());
        assertEquals(location2, actionableTrials.get(2).cancerType());
        assertEquals(loc2Doid1, actionableTrials.get(2).doid());
    }

    @Test
    public void canExtractAndMapDoid() {
        assertEquals("0060463", ActionableTrialFactory.extractDoid("0060463"));
        assertEquals("162", ActionableTrialFactory.extractDoid("0050686"));
    }
}