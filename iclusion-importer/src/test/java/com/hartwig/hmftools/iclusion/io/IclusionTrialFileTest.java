package com.hartwig.hmftools.iclusion.io;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.iclusion.data.IclusionMutation;
import com.hartwig.hmftools.iclusion.data.IclusionMutationLogicType;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionMutation;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionMutationCondition;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionTrial;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionTumorLocation;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class IclusionTrialFileTest {

    private static final IclusionMutation MUTATION1 =
            ImmutableIclusionMutation.builder().gene("gene1").name("name1").negation(false).build();

    private static final IclusionMutation MUTATION2 =
            ImmutableIclusionMutation.builder().gene("gene2").name("name2").negation(false).build();

    private static final IclusionMutation MUTATION3 =
            ImmutableIclusionMutation.builder().gene("gene3").name("name3").negation(false).build();

    @Test
    public void testEmptyTrialsConversion() {
        assertTrue(IclusionTrialFile.fromLines(IclusionTrialFile.toLines(Lists.newArrayList())).isEmpty());
    }

    @Test
    public void testRealTrialsConversionFromAndToString() {
        IclusionTrial trial1 = ImmutableIclusionTrial.builder()
                .id("1")
                .acronym("acr1")
                .title("title1")
                .eudra("eudra1")
                .nct("nct1")
                .ipn("ipn1")
                .ccmo("ccmo1")
                .addTumorLocations(ImmutableIclusionTumorLocation.builder()
                        .primaryTumorLocation("loc1")
                        .addDoids("doid1")
                        .addDoids("doid2")
                        .build())
                .addTumorLocations(ImmutableIclusionTumorLocation.builder().primaryTumorLocation("loc2").build())
                .addBlacklistedTumorLocations(ImmutableIclusionTumorLocation.builder().primaryTumorLocation("loc3").build())
                .addMutationConditions(ImmutableIclusionMutationCondition.builder()
                        .addMutations(MUTATION1)
                        .logicType(IclusionMutationLogicType.OR)
                        .build())
                .addMutationConditions(ImmutableIclusionMutationCondition.builder()
                        .addMutations(MUTATION2)
                        .logicType(IclusionMutationLogicType.OR)
                        .build())
                .build();

        IclusionTrial trial2 = ImmutableIclusionTrial.builder()
                .id("2")
                .acronym("acr2")
                .title("title2")
                .eudra("eudra2")
                .nct("nct2")
                .ipn("ipn2")
                .ccmo("ccmo2")
                .build();

        IclusionTrial trial3 = ImmutableIclusionTrial.builder()
                .id("3")
                .acronym("acr3")
                .title("title3")
                .eudra("eudra3")
                .nct("nct3")
                .ipn("ipn3")
                .ccmo("ccmo3")
                .addMutationConditions(ImmutableIclusionMutationCondition.builder()
                        .addMutations(MUTATION3)
                        .logicType(IclusionMutationLogicType.OR)
                        .build())
                .build();

        List<IclusionTrial> convertedTrials =
                IclusionTrialFile.fromLines(IclusionTrialFile.toLines(Lists.newArrayList(trial1, trial2, trial3)));

        assertEquals(3, convertedTrials.size());

        assertTrialsAreEqual(trial1, convertedTrials.get(0));
        assertTrialsAreEqual(trial2, convertedTrials.get(1));
        assertTrialsAreEqual(trial3, convertedTrials.get(2));
    }

    @Test
    public void ignoreLocationsAndMutationsWithInvalidCharacters() {
        IclusionTrial trialWithInvalidTumorLocation =
                testBuilder().addTumorLocations(ImmutableIclusionTumorLocation.builder().primaryTumorLocation("loc1 - loc2").build())
                        .build();

        IclusionMutation invalidMutation = ImmutableIclusionMutation.builder().gene("gene1 - 2").name("name3").negation(false).build();
        IclusionTrial trialWithInvalidMutation = testBuilder().addMutationConditions(ImmutableIclusionMutationCondition.builder()
                .addMutations(invalidMutation)
                .logicType(IclusionMutationLogicType.OR)
                .build()).build();

        List<IclusionTrial> convertedTrials = IclusionTrialFile.fromLines(IclusionTrialFile.toLines(Lists.newArrayList(
                trialWithInvalidTumorLocation,
                trialWithInvalidMutation)));

        assertEquals(2, convertedTrials.size());
        IclusionTrial expectedTrial = testBuilder().build();
        assertTrialsAreEqual(expectedTrial, convertedTrials.get(0));
        assertTrialsAreEqual(expectedTrial, convertedTrials.get(1));
    }

    private static void assertTrialsAreEqual(@NotNull IclusionTrial expectedTrial, @NotNull IclusionTrial actualTrial) {
        assertEquals(expectedTrial.id(), actualTrial.id());
        assertEquals(expectedTrial.acronym(), actualTrial.acronym());
        assertEquals(expectedTrial.title(), actualTrial.title());
        assertEquals(expectedTrial.eudra(), actualTrial.eudra());
        assertEquals(expectedTrial.nct(), actualTrial.nct());
        assertEquals(expectedTrial.ipn(), actualTrial.ipn());
        assertEquals(expectedTrial.ccmo(), actualTrial.ccmo());

        assertEquals(expectedTrial.tumorLocations().size(), actualTrial.tumorLocations().size());
        for (int i = 0; i < expectedTrial.tumorLocations().size(); i++) {
            assertEquals(expectedTrial.tumorLocations().get(i).primaryTumorLocation(),
                    actualTrial.tumorLocations().get(i).primaryTumorLocation());
            assertEquals(expectedTrial.tumorLocations().get(i).doids(), actualTrial.tumorLocations().get(i).doids());
        }

        assertEquals(expectedTrial.blacklistedTumorLocations().size(), actualTrial.blacklistedTumorLocations().size());
        for (int i = 0; i < expectedTrial.blacklistedTumorLocations().size(); i++) {
            assertEquals(expectedTrial.blacklistedTumorLocations().get(i).primaryTumorLocation(),
                    actualTrial.blacklistedTumorLocations().get(i).primaryTumorLocation());
            assertEquals(expectedTrial.blacklistedTumorLocations().get(i).doids(), actualTrial.blacklistedTumorLocations().get(i).doids());
        }

        assertEquals(expectedTrial.mutationConditions().size(), actualTrial.mutationConditions().size());
        for (int i = 0; i < expectedTrial.mutationConditions().size(); i++) {
            List<IclusionMutation> expectedMutations = expectedTrial.mutationConditions().get(i).mutations();
            List<IclusionMutation> actualMutations = actualTrial.mutationConditions().get(i).mutations();
            assertEquals(expectedMutations.size(), actualMutations.size());
            for (int j = 0; j < expectedMutations.size(); j++) {
                assertEquals(expectedMutations.get(j).gene(), actualMutations.get(j).gene());
                assertEquals(expectedMutations.get(j).name(), actualMutations.get(j).name());
                assertEquals(expectedMutations.get(j).negation(), actualMutations.get(j).negation());
            }
            assertEquals(expectedTrial.mutationConditions().get(i).logicType(), actualTrial.mutationConditions().get(i).logicType());
        }
    }

    @NotNull
    private static ImmutableIclusionTrial.Builder testBuilder() {
        return ImmutableIclusionTrial.builder()
                .id(Strings.EMPTY)
                .acronym(Strings.EMPTY)
                .title(Strings.EMPTY)
                .eudra(Strings.EMPTY)
                .nct(Strings.EMPTY)
                .ipn(Strings.EMPTY)
                .ccmo(Strings.EMPTY);
    }
}