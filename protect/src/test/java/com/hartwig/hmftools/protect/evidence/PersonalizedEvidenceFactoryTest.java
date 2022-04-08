package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.doid.DoidEdge;
import com.hartwig.hmftools.common.doid.DoidParents;
import com.hartwig.hmftools.common.doid.DoidParentsTest;
import com.hartwig.hmftools.common.protect.ProtectEvidenceType;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.ServeTestFactory;
import com.hartwig.hmftools.serve.actionability.ActionabilityTestUtil;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.actionability.characteristic.ImmutableActionableCharacteristic;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.hotspot.ImmutableActionableHotspot;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.actionability.range.ImmutableActionableRange;
import com.hartwig.hmftools.serve.actionability.range.RangeType;
import com.hartwig.hmftools.serve.cancertype.ImmutableCancerType;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicAnnotation;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicsComparator;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PersonalizedEvidenceFactoryTest {

    private static final double EPSILON = 1e-10;

    @NotNull
    public PersonalizedEvidenceFactory createPersonalizedEvidenceFactory(@NotNull String patientDoid) {

        Set<String> patientTumorDoids = Sets.newHashSet(patientDoid);
        return new PersonalizedEvidenceFactory(patientTumorDoids);
    }

    @NotNull
    public ActionableHotspot createActionableHotspot(@NotNull String blacklistCancerType, @NotNull String blacklistDoid,
            @NotNull String cancerType, @NotNull String doid) {
        ActionableEvent event = ActionabilityTestUtil.create(Knowledgebase.CKB,
                "amp",
                Sets.newHashSet(),
                "treatment A",
                ImmutableCancerType.builder().name(cancerType).doid(doid).build(),
                Sets.newHashSet(ImmutableCancerType.builder().name(blacklistCancerType).doid(blacklistDoid).build()),
                EvidenceLevel.A,
                EvidenceDirection.RESPONSIVE,
                Sets.newHashSet());

        return ImmutableActionableHotspot.builder()
                .from(event)
                .chromosome(Strings.EMPTY)
                .position(0)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .build();
    }

    @Test
    public void canDetermineOnlabel() {
        PersonalizedEvidenceFactory factoryOnlabelMatching = createPersonalizedEvidenceFactory("162");
        ActionableHotspot hotspotOnLabelMatching = createActionableHotspot("", "", "Cancer", "162");
        assertTrue(factoryOnlabelMatching.determineOnlabel(hotspotOnLabelMatching.applicableCancerType(),
                hotspotOnLabelMatching.blacklistCancerTypes()));

        PersonalizedEvidenceFactory factoryNotBlacklisted = createPersonalizedEvidenceFactory("10283");
        ActionableHotspot hotspotNotBlacklisted = createActionableHotspot("Breast", "0060081", "prostate", "10283");
        assertTrue(factoryNotBlacklisted.determineOnlabel(hotspotNotBlacklisted.applicableCancerType(),
                hotspotNotBlacklisted.blacklistCancerTypes()));

        PersonalizedEvidenceFactory factoryBlacklisted = createPersonalizedEvidenceFactory("10283");
        ActionableHotspot hotspotBlacklisted = createActionableHotspot("Prostate", "10283", "Cancer", "162");
        assertFalse(factoryBlacklisted.determineOnlabel(hotspotBlacklisted.applicableCancerType(),
                hotspotBlacklisted.blacklistCancerTypes()));
    }

    @Test
    public void canDetermineBlacklistedEvidence() {
        PersonalizedEvidenceFactory factoryBlacklisted = createPersonalizedEvidenceFactory("10283");
        ActionableHotspot hotspotBlacklisted = createActionableHotspot("Prostate", "10283", "Cancer", "162");
        assertTrue(factoryBlacklisted.determineBlacklistedEvidence(hotspotBlacklisted.blacklistCancerTypes()));

        PersonalizedEvidenceFactory factoryNotMatchWithBlacklisted = createPersonalizedEvidenceFactory("0060081");
        ActionableHotspot hotspotNotMatchWithBlacklisted = createActionableHotspot("Prostate", "10283", "Cancer", "162");
        assertFalse(factoryNotMatchWithBlacklisted.determineBlacklistedEvidence(hotspotNotMatchWithBlacklisted.blacklistCancerTypes()));

        PersonalizedEvidenceFactory factoryNotBlacklisted = createPersonalizedEvidenceFactory("10383");
        ActionableHotspot hotspotNotBlacklisted = createActionableHotspot("", "", "Cancer", "162");
        assertFalse(factoryNotBlacklisted.determineBlacklistedEvidence(hotspotNotBlacklisted.blacklistCancerTypes()));
    }

    @Test
    public void canDetermineEvidenceTypes() {
        assertEquals(ProtectEvidenceType.HOTSPOT_MUTATION,
                PersonalizedEvidenceFactory.determineEvidenceType(ServeTestFactory.createTestActionableHotspot()));

        ActionableRange range =
                ImmutableActionableRange.builder().from(ServeTestFactory.createTestActionableRange()).rangeType(RangeType.EXON).build();
        assertEquals(ProtectEvidenceType.EXON_MUTATION, PersonalizedEvidenceFactory.determineEvidenceType(range));

        ActionableGene gene = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .event(GeneLevelEvent.INACTIVATION)
                .build();
        assertEquals(ProtectEvidenceType.INACTIVATION, PersonalizedEvidenceFactory.determineEvidenceType(gene));

        assertEquals(ProtectEvidenceType.FUSION_PAIR,
                PersonalizedEvidenceFactory.determineEvidenceType(ServeTestFactory.createTestActionableFusion()));

        ActionableCharacteristic characteristic = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic(null, null))
                .name(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD)
                .comparator(TumorCharacteristicsComparator.EQUAL_OR_GREATER)
                .cutoff((double) 10)
                .build();
        assertEquals(ProtectEvidenceType.SIGNATURE, PersonalizedEvidenceFactory.determineEvidenceType(characteristic));
        assertEquals(10, characteristic.cutoff(), EPSILON);
    }

    @Test
    public void canDetermineEvidenceTypesForAllRanges() {
        for (RangeType rangeType : RangeType.values()) {
            ActionableRange range =
                    ImmutableActionableRange.builder().from(ServeTestFactory.createTestActionableRange()).rangeType(rangeType).build();
            assertNotNull(PersonalizedEvidenceFactory.determineEvidenceType(range));
        }
    }

    @Test
    public void canDetermineEvidenceTypesFroAllGeneEvents() {
        for (GeneLevelEvent geneLevelEvent : GeneLevelEvent.values()) {
            ActionableGene gene =
                    ImmutableActionableGene.builder().from(ServeTestFactory.createTestActionableGene()).event(geneLevelEvent).build();
            assertNotNull(PersonalizedEvidenceFactory.determineEvidenceType(gene));
        }
    }

    @Test
    public void canDetermineEvidenceTypesForAllCharacteristics() {
        for (TumorCharacteristicAnnotation name : TumorCharacteristicAnnotation.values()) {
            ActionableCharacteristic characteristic = ImmutableActionableCharacteristic.builder()
                    .from(ServeTestFactory.createTestActionableCharacteristic(null, null))
                    .name(name)
                    .build();
            assertNotNull(PersonalizedEvidenceFactory.determineEvidenceType(characteristic));
        }
    }

    @Test
    public void canDetermineRangeRank() {
        ActionableRange range = ImmutableActionableRange.builder().from(ServeTestFactory.createTestActionableRange()).rank(2).build();

        assertEquals(2, (int) PersonalizedEvidenceFactory.determineRangeRank(range));

        assertNull(PersonalizedEvidenceFactory.determineRangeRank(ServeTestFactory.createTestActionableFusion()));
    }
}