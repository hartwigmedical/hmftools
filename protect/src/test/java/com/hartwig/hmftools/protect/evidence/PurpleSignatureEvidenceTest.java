package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.ImmutablePurpleData;
import com.hartwig.hmftools.common.purple.PurpleData;
import com.hartwig.hmftools.common.purple.PurpleTestFactory;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.serve.ServeTestFactory;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.actionability.characteristic.ImmutableActionableCharacteristic;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicAnnotation;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicsComparator;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleSignatureEvidenceTest {

    @Test
    public void canHandleNonPurpleSignatureEvidence() {
        ActionableCharacteristic nonPurple = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic())
                .name(TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT)
                .build();

        PurpleSignatureEvidence purpleSignatureEvidence = new PurpleSignatureEvidence(EvidenceTestFactory.createTestEvidenceFactory(),
                Lists.newArrayList(nonPurple));

        PurpleData data = createPurpleData(MicrosatelliteStatus.MSI, TumorMutationalStatus.HIGH, 4.0, 241, 0.0);
        assertEquals(0, purpleSignatureEvidence.evidence(data).size());
    }

    @Test
    public void canDetermineMSI() {
        ActionableCharacteristic signatureDefault = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic())
                .name(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE)
                .build();

        ActionableCharacteristic signatureWithCutoff = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic())
                .name(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE)
                .comparator(TumorCharacteristicsComparator.EQUAL_OR_GREATER)
                .cutoff(4D)
                .build();

        PurpleSignatureEvidence purpleSignatureEvidence = new PurpleSignatureEvidence(EvidenceTestFactory.createTestEvidenceFactory(),
                Lists.newArrayList(signatureDefault, signatureWithCutoff));

        PurpleData both = createPurpleData(MicrosatelliteStatus.MSI, TumorMutationalStatus.HIGH, 4.0, 241, 0.0);
        assertEquals(2, purpleSignatureEvidence.evidence(both).size());

        PurpleData one = createPurpleData(MicrosatelliteStatus.MSI, TumorMutationalStatus.HIGH, 3.0, 241, 0.0);
        assertEquals(1, purpleSignatureEvidence.evidence(one).size());

        PurpleData none = createPurpleData(MicrosatelliteStatus.MSS, TumorMutationalStatus.HIGH, 3.0, 241, 0.0);
        assertEquals(0, purpleSignatureEvidence.evidence(none).size());
    }

    @Test
    public void canDetermineMSS() {
        ActionableCharacteristic signatureDefault = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic())
                .name(TumorCharacteristicAnnotation.MICROSATELLITE_STABLE)
                .build();

        ActionableCharacteristic signatureWithCutoff = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic())
                .name(TumorCharacteristicAnnotation.MICROSATELLITE_STABLE)
                .comparator(TumorCharacteristicsComparator.LOWER)
                .cutoff(4D)
                .build();

        PurpleSignatureEvidence purpleSignatureEvidence = new PurpleSignatureEvidence(EvidenceTestFactory.createTestEvidenceFactory(),
                Lists.newArrayList(signatureDefault, signatureWithCutoff));

        PurpleData both = createPurpleData(MicrosatelliteStatus.MSS, TumorMutationalStatus.HIGH, 2.0, 241, 0.0);
        assertEquals(2, purpleSignatureEvidence.evidence(both).size());

        PurpleData one = createPurpleData(MicrosatelliteStatus.MSS, TumorMutationalStatus.HIGH, 5.0, 241, 0.0);
        assertEquals(1, purpleSignatureEvidence.evidence(one).size());

        PurpleData none = createPurpleData(MicrosatelliteStatus.MSI, TumorMutationalStatus.HIGH, 5.0, 241, 0.0);
        assertEquals(0, purpleSignatureEvidence.evidence(none).size());
    }

    @Test
    public void canDetermineHighTML() {
        ActionableCharacteristic signatureDefault = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic())
                .name(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD)
                .build();

        ActionableCharacteristic signatureWithCutoff = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic())
                .name(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD)
                .comparator(TumorCharacteristicsComparator.GREATER)
                .cutoff(100D)
                .build();

        PurpleSignatureEvidence purpleSignatureEvidence = new PurpleSignatureEvidence(EvidenceTestFactory.createTestEvidenceFactory(),
                Lists.newArrayList(signatureDefault, signatureWithCutoff));

        PurpleData both = createPurpleData(MicrosatelliteStatus.MSS, TumorMutationalStatus.HIGH, 2.0, 241, 0.0);
        assertEquals(2, purpleSignatureEvidence.evidence(both).size());

        PurpleData one = createPurpleData(MicrosatelliteStatus.MSS, TumorMutationalStatus.LOW, 5.0, 120, 0.0);
        assertEquals(1, purpleSignatureEvidence.evidence(one).size());

        PurpleData none = createPurpleData(MicrosatelliteStatus.MSI, TumorMutationalStatus.LOW, 5.0, 80, 0.0);
        assertEquals(0, purpleSignatureEvidence.evidence(none).size());
    }

    @Test
    public void canDetermineLowTML() {
        ActionableCharacteristic signatureDefault = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic())
                .name(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD)
                .build();

        ActionableCharacteristic signatureWithCutoff = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic())
                .name(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD)
                .comparator(TumorCharacteristicsComparator.EQUAL_OR_LOWER)
                .cutoff(100D)
                .build();

        PurpleSignatureEvidence purpleSignatureEvidence = new PurpleSignatureEvidence(EvidenceTestFactory.createTestEvidenceFactory(),
                Lists.newArrayList(signatureDefault, signatureWithCutoff));

        PurpleData both = createPurpleData(MicrosatelliteStatus.MSS, TumorMutationalStatus.LOW, 2.0, 80, 0.0);
        assertEquals(2, purpleSignatureEvidence.evidence(both).size());

        PurpleData one = createPurpleData(MicrosatelliteStatus.MSS, TumorMutationalStatus.HIGH, 5.0, 100, 0.0);
        assertEquals(1, purpleSignatureEvidence.evidence(one).size());

        PurpleData none = createPurpleData(MicrosatelliteStatus.MSI, TumorMutationalStatus.HIGH, 5.0, 120, 0.0);
        assertEquals(0, purpleSignatureEvidence.evidence(none).size());
    }

    @Test
    public void canDetermineHighTMB() {
        ActionableCharacteristic signatureDefault = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic())
                .name(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_BURDEN)
                .build();

        ActionableCharacteristic signatureWithCutoff = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic())
                .name(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_BURDEN)
                .comparator(TumorCharacteristicsComparator.EQUAL_OR_GREATER)
                .cutoff(12D)
                .build();

        PurpleSignatureEvidence purpleSignatureEvidence = new PurpleSignatureEvidence(EvidenceTestFactory.createTestEvidenceFactory(),
                Lists.newArrayList(signatureDefault, signatureWithCutoff));

        PurpleData both = createPurpleData(MicrosatelliteStatus.MSS, TumorMutationalStatus.LOW, 2.0, 80, 12.0);
        assertEquals(2, purpleSignatureEvidence.evidence(both).size());

        PurpleData one = createPurpleData(MicrosatelliteStatus.MSS, TumorMutationalStatus.HIGH, 5.0, 100, 11.0);
        assertEquals(1, purpleSignatureEvidence.evidence(one).size());

        PurpleData none = createPurpleData(MicrosatelliteStatus.MSI, TumorMutationalStatus.HIGH, 5.0, 120, 8.0);
        assertEquals(0, purpleSignatureEvidence.evidence(none).size());
    }

    @Test
    public void canDetermineLowTMB() {
        ActionableCharacteristic signatureDefault = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic())
                .name(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_BURDEN)
                .build();

        ActionableCharacteristic signatureWithCutoff = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic())
                .name(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_BURDEN)
                .comparator(TumorCharacteristicsComparator.LOWER)
                .cutoff(8D)
                .build();

        PurpleSignatureEvidence purpleSignatureEvidence = new PurpleSignatureEvidence(EvidenceTestFactory.createTestEvidenceFactory(),
                Lists.newArrayList(signatureDefault, signatureWithCutoff));

        PurpleData both = createPurpleData(MicrosatelliteStatus.MSS, TumorMutationalStatus.LOW, 2.0, 80, 7D);
        assertEquals(2, purpleSignatureEvidence.evidence(both).size());

        PurpleData one = createPurpleData(MicrosatelliteStatus.MSS, TumorMutationalStatus.HIGH, 5.0, 100, 9D);
        assertEquals(1, purpleSignatureEvidence.evidence(one).size());

        PurpleData none = createPurpleData(MicrosatelliteStatus.MSI, TumorMutationalStatus.HIGH, 5.0, 120, 11D);
        assertEquals(0, purpleSignatureEvidence.evidence(none).size());
    }


    @Test
    public void canConvertCharacteristicToEvent() {
        assertEquals("Microsatellite unstable", PurpleSignatureEvidence.toEvent(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE));
    }

    @NotNull
    private static PurpleData createPurpleData(@NotNull MicrosatelliteStatus msStatus, @NotNull TumorMutationalStatus tmlStatus,
            double microsatelliteIndelsPerMb, int tumorMutationalLoad, double tumorMutationalBurdenPerMb) {
        return ImmutablePurpleData.builder()
                .from(PurpleTestFactory.createMinimalTestPurpleData())
                .microsatelliteStatus(msStatus)
                .microsatelliteIndelsPerMb(microsatelliteIndelsPerMb)
                .tumorMutationalLoadStatus(tmlStatus)
                .tumorMutationalLoad(tumorMutationalLoad)
                .tumorMutationalBurdenPerMb(tumorMutationalBurdenPerMb)
                .build();
    }
}