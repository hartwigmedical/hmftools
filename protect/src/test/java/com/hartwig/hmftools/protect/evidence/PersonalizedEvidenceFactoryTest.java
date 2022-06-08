package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.google.common.collect.Sets;
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
import com.hartwig.hmftools.serve.cancertype.CancerType;
import com.hartwig.hmftools.serve.cancertype.ImmutableCancerType;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicAnnotation;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.treatment.ImmutableTreatment;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class PersonalizedEvidenceFactoryTest {

    @Test
    public void canDetermineOnLabel() {
        PersonalizedEvidenceFactory factoryOnLabelMatching = EvidenceTestFactory.create("162");
        ActionableHotspot hotspotOnLabelMatching = create("Cancer", "162");
        assertTrue(factoryOnLabelMatching.isOnLabel(hotspotOnLabelMatching.applicableCancerType(),
                hotspotOnLabelMatching.blacklistCancerTypes(),
                "treatment"));

        PersonalizedEvidenceFactory factoryNotBlacklisted = EvidenceTestFactory.create("10283");
        ActionableHotspot hotspotNotBlacklisted = create("prostate", "10283", "Breast", "0060081");
        assertTrue(factoryNotBlacklisted.isOnLabel(hotspotNotBlacklisted.applicableCancerType(),
                hotspotNotBlacklisted.blacklistCancerTypes(),
                "treatment"));

        PersonalizedEvidenceFactory factoryBlacklisted = EvidenceTestFactory.create("10283");
        ActionableHotspot hotspotBlacklisted = create("Cancer", "162", "Prostate", "10283");
        assertFalse(factoryBlacklisted.isOnLabel(hotspotBlacklisted.applicableCancerType(),
                hotspotBlacklisted.blacklistCancerTypes(),
                "treatment"));
    }

    @Test
    public void canDetermineBlacklistedEvidence() {
        PersonalizedEvidenceFactory factoryBlacklisted = EvidenceTestFactory.create("10283");
        ActionableHotspot hotspotBlacklisted = create("Cancer", "162", "Prostate", "10283");
        assertTrue(factoryBlacklisted.determineBlacklistedEvidence(hotspotBlacklisted.blacklistCancerTypes(), "treatment"));

        PersonalizedEvidenceFactory factoryNotMatchWithBlacklisted = EvidenceTestFactory.create("0060081");
        ActionableHotspot hotspotNotMatchWithBlacklisted = create("Cancer", "162", "Prostate", "10283");
        assertFalse(factoryNotMatchWithBlacklisted.determineBlacklistedEvidence(hotspotNotMatchWithBlacklisted.blacklistCancerTypes(),
                "treatment"));

        PersonalizedEvidenceFactory factoryNotBlacklisted = EvidenceTestFactory.create("10383");
        ActionableHotspot hotspotNotBlacklisted = create("Cancer", "162");
        assertFalse(factoryNotBlacklisted.determineBlacklistedEvidence(hotspotNotBlacklisted.blacklistCancerTypes(), "treatment"));
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

        ActionableGene amplification = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .event(GeneLevelEvent.AMPLIFICATION)
                .build();
        assertEquals(ProtectEvidenceType.AMPLIFICATION, PersonalizedEvidenceFactory.determineEvidenceType(amplification));

        ActionableGene overexpression = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .event(GeneLevelEvent.OVER_EXPRESSION)
                .build();
        assertEquals(ProtectEvidenceType.OVER_EXPRESSION, PersonalizedEvidenceFactory.determineEvidenceType(overexpression));


        assertEquals(ProtectEvidenceType.FUSION_PAIR,
                PersonalizedEvidenceFactory.determineEvidenceType(ServeTestFactory.createTestActionableFusion()));

        assertEquals(ProtectEvidenceType.SIGNATURE,
                PersonalizedEvidenceFactory.determineEvidenceType(ServeTestFactory.createTestActionableCharacteristic()));
    }

    @Test
    public void canDetermineEvidenceTypesForAllRanges() {
        ActionableRange base = ServeTestFactory.createTestActionableRange();
        for (RangeType rangeType : RangeType.values()) {
            ActionableRange range = ImmutableActionableRange.builder().from(base).rangeType(rangeType).build();
            assertNotNull(PersonalizedEvidenceFactory.determineEvidenceType(range));
        }
    }

    @Test
    public void canDetermineEvidenceTypesFroAllGeneEvents() {
        ActionableGene base = ServeTestFactory.createTestActionableGene();
        for (GeneLevelEvent geneLevelEvent : GeneLevelEvent.values()) {
            ActionableGene gene = ImmutableActionableGene.builder().from(base).event(geneLevelEvent).build();
            assertNotNull(PersonalizedEvidenceFactory.determineEvidenceType(gene));
        }
    }

    @Test
    public void canDetermineEvidenceTypesForAllCharacteristics() {
        ActionableCharacteristic base = ServeTestFactory.createTestActionableCharacteristic();
        for (TumorCharacteristicAnnotation name : TumorCharacteristicAnnotation.values()) {
            ActionableCharacteristic characteristic = ImmutableActionableCharacteristic.builder().from(base).name(name).build();
            assertNotNull(PersonalizedEvidenceFactory.determineEvidenceType(characteristic));
        }
    }

    @Test
    public void canDetermineRangeRank() {
        ActionableRange range = ImmutableActionableRange.builder().from(ServeTestFactory.createTestActionableRange()).rank(2).build();

        assertEquals(2, (int) PersonalizedEvidenceFactory.determineRangeRank(range));

        assertNull(PersonalizedEvidenceFactory.determineRangeRank(ServeTestFactory.createTestActionableFusion()));
    }

    @NotNull
    private static ActionableHotspot create(@NotNull String cancerType, @NotNull String doid) {
        return create(cancerType, doid, null, null);
    }

    @NotNull
    private static ActionableHotspot create(@NotNull String cancerType, @NotNull String doid, @Nullable String blacklistCancerType,
            @Nullable String blacklistDoid) {
        Set<CancerType> blacklist = Sets.newHashSet();
        if (blacklistCancerType != null && blacklistDoid != null) {
            blacklist.add(ImmutableCancerType.builder().name(blacklistCancerType).doid(blacklistDoid).build());
        }

        ActionableEvent event = ActionabilityTestUtil.create(Knowledgebase.CKB,
                "amp",
                Sets.newHashSet(),
                ImmutableTreatment.builder().treament("treatment A").drugClasses(Sets.newHashSet("drugClasses")).build(),
                ImmutableCancerType.builder().name(cancerType).doid(doid).build(),
                blacklist,
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
}