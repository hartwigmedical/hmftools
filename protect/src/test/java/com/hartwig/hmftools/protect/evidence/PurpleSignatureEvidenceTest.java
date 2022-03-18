package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
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
    public void canDeterminePurpleSignatureEvidence() {
        ActionableCharacteristic signature1 = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic(TumorCharacteristicsComparator.GREATER, 240.0))
                .name(TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT)
                .build();

        ActionableCharacteristic signature2 = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic(TumorCharacteristicsComparator.GREATER, 240.0))
                .name(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD)
                .build();

        ActionableCharacteristic signature3 = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic(TumorCharacteristicsComparator.EQUAL_OR_GREATER, 4.0))
                .name(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE)
                .build();

        ActionableCharacteristic signature4 = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic(TumorCharacteristicsComparator.EQUAL_OR_GREATER, 11.0))
                .name(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_BURDEN)
                .build();

        PurpleSignatureEvidence purpleSignatureEvidence = new PurpleSignatureEvidence(EvidenceTestFactory.createTestEvidenceFactory(),
                Lists.newArrayList(signature1, signature2, signature3, signature4));

        PurpleData all = createPurpleData(MicrosatelliteStatus.MSI, TumorMutationalStatus.HIGH, 4.0, 241, 0.0);
        List<ProtectEvidence> evidences = purpleSignatureEvidence.evidence(all);
               assertEquals(2, evidences.size());
        ProtectEvidence msiEvidence = find(evidences, PurpleSignatureEvidence.MICROSATELLITE_UNSTABLE_EVENT);
        assertTrue(msiEvidence.reported());

        ProtectEvidence tmlEvidence = find(evidences, PurpleSignatureEvidence.HIGH_TUMOR_LOAD_EVENT);
        assertTrue(tmlEvidence.reported());

        PurpleData msi = createPurpleData(MicrosatelliteStatus.MSI, TumorMutationalStatus.LOW, 4.0, 0, 0.0);
        assertEquals(1, purpleSignatureEvidence.evidence(msi).size());

        PurpleData tmlHigh = createPurpleData(MicrosatelliteStatus.MSS, TumorMutationalStatus.HIGH, 0.0, 241, 0.0);
        assertEquals(1, purpleSignatureEvidence.evidence(tmlHigh).size());

        PurpleData tmbHigh = createPurpleData(MicrosatelliteStatus.MSS, TumorMutationalStatus.HIGH, 0.0, 240, 11.0);
        assertEquals(1, purpleSignatureEvidence.evidence(tmbHigh).size());

        PurpleData none = createPurpleData(MicrosatelliteStatus.MSS, TumorMutationalStatus.LOW, 0.0, 240, 0.0);
        assertTrue(purpleSignatureEvidence.evidence(none).isEmpty());
    }

    @NotNull
    private static ProtectEvidence find(@NotNull List<ProtectEvidence> evidences, @NotNull String eventToFind) {
        for (ProtectEvidence evidence : evidences) {
            if (evidence.event().equals(eventToFind)) {
                return evidence;
            }
        }

        throw new IllegalStateException("Could not event in evidences: " + eventToFind);
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