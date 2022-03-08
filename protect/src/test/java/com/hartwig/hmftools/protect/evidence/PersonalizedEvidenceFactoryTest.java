package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.common.protect.ProtectEvidenceType;
import com.hartwig.hmftools.serve.ServeTestFactory;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.actionability.characteristic.ImmutableActionableCharacteristic;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.actionability.range.ImmutableActionableRange;
import com.hartwig.hmftools.serve.actionability.range.RangeType;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicAnnotation;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicsComparator;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;

import org.junit.Test;

public class PersonalizedEvidenceFactoryTest {

    private static final double EPSILON = 1e-10;

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
                .from(ServeTestFactory.createTestActionableCharacteristic())
                .name(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD)
                .comparator(TumorCharacteristicsComparator.EQUAL_OR_GREATER)
                .cutoff((double)10)
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
                    .from(ServeTestFactory.createTestActionableCharacteristic())
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