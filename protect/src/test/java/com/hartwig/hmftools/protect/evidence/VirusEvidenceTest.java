package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.ImmutableVirusInterpreterData;
import com.hartwig.hmftools.common.virus.VirusConstants;
import com.hartwig.hmftools.common.virus.VirusInterpreterData;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;
import com.hartwig.hmftools.common.virus.VirusTestFactory;
import com.hartwig.hmftools.serve.ServeTestFactory;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.actionability.characteristic.ImmutableActionableCharacteristic;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicAnnotation;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class VirusEvidenceTest {

    @Test
    public void canDetermineEvidenceForViruses() {
        VirusInterpreterData testData = createTestVirusInterpreterData();

        ActionableCharacteristic hpv = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic(null, null))
                .name(TumorCharacteristicAnnotation.HPV_POSITIVE)
                .build();

        ActionableCharacteristic ebv = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic(null, null))
                .name(TumorCharacteristicAnnotation.EBV_POSITIVE)
                .build();

        VirusEvidence virusEvidence = new VirusEvidence(EvidenceTestFactory.createTestEvidenceFactory(), Lists.newArrayList(hpv, ebv));

        List<ProtectEvidence> evidences = virusEvidence.evidence(testData);
        assertEquals(2, evidences.size());

        // The test data has a reportable HPV virus
        ProtectEvidence hpvEvidence = find(evidences, VirusEvidence.HPV_POSITIVE_EVENT);
        assertTrue(hpvEvidence.reported());

        // The test data has a reportable LOW driver EBV virus
        ProtectEvidence ebvEvidence = find(evidences, VirusEvidence.EBV_POSITIVE_EVENT);
        assertFalse(ebvEvidence.reported());
    }

    @NotNull
    private static ProtectEvidence find(@NotNull List<ProtectEvidence> evidences, @NotNull String event) {
        for (ProtectEvidence evidence : evidences) {
            if (evidence.event().equals(event)) {
                return evidence;
            }
        }

        throw new IllegalStateException("Could not find evidence with genomic event: " + event);
    }

    @NotNull
    private static VirusInterpreterData createTestVirusInterpreterData() {
        List<AnnotatedVirus> reportable = Lists.newArrayList();
        reportable.add(VirusTestFactory.testAnnotatedVirusBuilder()
                .interpretation(VirusConstants.HPV.name())
                .reported(true)
                .virusDriverLikelihoodType(VirusLikelihoodType.HIGH)
                .build());
        reportable.add(VirusTestFactory.testAnnotatedVirusBuilder()
                .interpretation(VirusConstants.MCV.name())
                .reported(true)
                .virusDriverLikelihoodType(VirusLikelihoodType.LOW)
                .build());
        reportable.add(VirusTestFactory.testAnnotatedVirusBuilder()
                .interpretation(VirusConstants.EBV.name())
                .reported(true)
                .virusDriverLikelihoodType(VirusLikelihoodType.LOW)
                .build());
        reportable.add(VirusTestFactory.testAnnotatedVirusBuilder()
                .reported(true)
                .virusDriverLikelihoodType(VirusLikelihoodType.UNKNOWN)
                .build());

        List<AnnotatedVirus> unreported = Lists.newArrayList();
        unreported.add(VirusTestFactory.testAnnotatedVirusBuilder()
                .interpretation(VirusConstants.EBV.name())
                .reported(false)
                .virusDriverLikelihoodType(VirusLikelihoodType.HIGH)
                .build());
        unreported.add(VirusTestFactory.testAnnotatedVirusBuilder()
                .interpretation(VirusConstants.EBV.name())
                .reported(false)
                .virusDriverLikelihoodType(VirusLikelihoodType.HIGH)
                .build());
        unreported.add(VirusTestFactory.testAnnotatedVirusBuilder()
                .reported(false)
                .virusDriverLikelihoodType(VirusLikelihoodType.UNKNOWN)
                .build());

        return ImmutableVirusInterpreterData.builder().reportableViruses(reportable).unreportedViruses(unreported).build();
    }
}