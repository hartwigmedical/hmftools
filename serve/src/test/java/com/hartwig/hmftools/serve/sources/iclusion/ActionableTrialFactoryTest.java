package com.hartwig.hmftools.serve.sources.iclusion;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTumorLocation;
import com.hartwig.hmftools.iclusion.datamodel.ImmutableIclusionTumorLocation;
import com.hartwig.hmftools.serve.cancertype.ImmutableCancerType;
import com.hartwig.hmftools.serve.curation.DoidLookupTestFactory;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class ActionableTrialFactoryTest {

    @Test
    public void canCreateActionableTrials() {
        String location1 = "loc1";
        String loc1Doid1 = "162";
        String loc1Doid2 = "loc1Doid2";
        String location2 = "loc2";
        String loc2Doid1 = "loc2Doid2";
        String blacklistLocation1 = "blacklistLocation";
        String blacklistDoid1 = "blacklistDoid";
        String treatment = "trial";

        IclusionTumorLocation loc1 =
                ImmutableIclusionTumorLocation.builder().primaryTumorLocation(location1).addDoids(loc1Doid1).addDoids(loc1Doid2).build();
        IclusionTumorLocation loc2 = ImmutableIclusionTumorLocation.builder().primaryTumorLocation(location2).addDoids(loc2Doid1).build();
        IclusionTumorLocation blacklist =
                ImmutableIclusionTumorLocation.builder().primaryTumorLocation(blacklistLocation1).addDoids(blacklistDoid1).build();
        IclusionTrial trial = IclusionTestFactory.trialWithTumors(treatment, Lists.newArrayList(loc1, loc2), Lists.newArrayList(blacklist));

        ActionableTrialFactory factory = new ActionableTrialFactory(DoidLookupTestFactory.dummy());
        List<ActionableTrial> actionableTrials = factory.toActionableTrials(trial, Strings.EMPTY);
        assertEquals(3, actionableTrials.size());
        assertEquals(treatment, actionableTrials.get(0).treatment());
        assertEquals(location1, actionableTrials.get(0).applicableCancerType().name());
        assertEquals(loc1Doid1, actionableTrials.get(0).applicableCancerType().doid());
        assertEquals(Sets.newHashSet(ImmutableCancerType.builder().name("blacklistLocation").doid("blacklistDoid").build(),
                        ImmutableCancerType.builder().name("Hematologic cancer").doid("2531").build()),
                actionableTrials.get(0).blacklistCancerTypes());

        assertEquals(treatment, actionableTrials.get(1).treatment());
        assertEquals(location1, actionableTrials.get(1).applicableCancerType().name());
        assertEquals(loc1Doid2, actionableTrials.get(1).applicableCancerType().doid());
        assertEquals(Sets.newHashSet(ImmutableCancerType.builder().name("blacklistLocation").doid("blacklistDoid").build(),
                        ImmutableCancerType.builder().name("Hematologic cancer").doid("2531").build()),
                actionableTrials.get(1).blacklistCancerTypes());

        assertEquals(treatment, actionableTrials.get(2).treatment());
        assertEquals(location2, actionableTrials.get(2).applicableCancerType().name());
        assertEquals(loc2Doid1, actionableTrials.get(2).applicableCancerType().doid());
        assertEquals(Sets.newHashSet(ImmutableCancerType.builder().name("blacklistLocation").doid("blacklistDoid").build(),
                        ImmutableCancerType.builder().name("Hematologic cancer").doid("2531").build()),
                actionableTrials.get(2).blacklistCancerTypes());

        IclusionTrial trialOnlyBlacklist =
                IclusionTestFactory.trialWithTumors(treatment, Lists.newArrayList(loc2), Lists.newArrayList(blacklist));
        List<ActionableTrial> actionableTrialsOnlyBlacklist = factory.toActionableTrials(trialOnlyBlacklist, Strings.EMPTY);
        assertEquals(1, actionableTrialsOnlyBlacklist.size());
        assertEquals(treatment, actionableTrialsOnlyBlacklist.get(0).treatment());
        assertEquals(location2, actionableTrialsOnlyBlacklist.get(0).applicableCancerType().name());
        assertEquals(loc2Doid1, actionableTrialsOnlyBlacklist.get(0).applicableCancerType().doid());
        assertEquals(Sets.newHashSet(ImmutableCancerType.builder().name("blacklistLocation").doid("blacklistDoid").build()),
                actionableTrialsOnlyBlacklist.get(0).blacklistCancerTypes());

        IclusionTrial trialWithoutBlacklist =
                IclusionTestFactory.trialWithTumors(treatment, Lists.newArrayList(loc2), Lists.newArrayList());
        List<ActionableTrial> actionableTrialsWithoutBlacklist = factory.toActionableTrials(trialWithoutBlacklist, Strings.EMPTY);
        assertEquals(1, actionableTrialsWithoutBlacklist.size());
        assertEquals(treatment, actionableTrialsWithoutBlacklist.get(0).treatment());
        assertEquals(location2, actionableTrialsWithoutBlacklist.get(0).applicableCancerType().name());
        assertEquals(loc2Doid1, actionableTrialsWithoutBlacklist.get(0).applicableCancerType().doid());
        assertEquals(Sets.newHashSet(), actionableTrialsWithoutBlacklist.get(0).blacklistCancerTypes());

        IclusionTrial trialWith162 = IclusionTestFactory.trialWithTumors(treatment, Lists.newArrayList(loc1), Lists.newArrayList());
        List<ActionableTrial> actionableTrialsWith162 = factory.toActionableTrials(trialWith162, Strings.EMPTY);
        assertEquals(2, actionableTrialsWith162.size());
        assertEquals(treatment, actionableTrialsWith162.get(0).treatment());
        assertEquals(location1, actionableTrialsWith162.get(0).applicableCancerType().name());
        assertEquals(loc1Doid1, actionableTrialsWith162.get(0).applicableCancerType().doid());
        assertEquals(Sets.newHashSet(ImmutableCancerType.builder().name("Hematologic cancer").doid("2531").build()),
                actionableTrialsWith162.get(0).blacklistCancerTypes());
    }

    @Test
    public void canExtractAndMapDoid() {
        assertEquals("0060463", ActionableTrialFactory.extractDoid("0060463"));
        assertEquals("162", ActionableTrialFactory.extractDoid("0050686"));
        assertEquals("162", ActionableTrialFactory.extractDoid("UNKNOWN"));
        assertEquals("162", ActionableTrialFactory.extractDoid("MESH: D009382"));
    }
}