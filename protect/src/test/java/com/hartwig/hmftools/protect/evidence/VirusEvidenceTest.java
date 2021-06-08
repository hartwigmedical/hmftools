package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusInterpretation;
import com.hartwig.hmftools.common.virus.VirusTestFactory;
import com.hartwig.hmftools.protect.ProtectTestFactory;
import com.hartwig.hmftools.protect.virusinterpreter.ImmutableVirusInterpreterData;
import com.hartwig.hmftools.protect.virusinterpreter.VirusInterpreterData;
import com.hartwig.hmftools.serve.ServeTestFactory;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.actionability.characteristic.ImmutableActionableCharacteristic;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristic;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class VirusEvidenceTest {

    @Test
    public void canDetermineEvidenceForViruses() {
        VirusInterpreterData testData = createTestVirusInterpreterData();

        ActionableCharacteristic hpvEvidence = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic())
                .name(TumorCharacteristic.HPV_POSITIVE)
                .build();

        ActionableCharacteristic ebvEvidence = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic())
                .name(TumorCharacteristic.EBV_POSITIVE)
                .build();

        VirusEvidence virusEvidence =
                new VirusEvidence(ProtectTestFactory.createTestEvidenceFactory(), Lists.newArrayList(hpvEvidence, ebvEvidence));

        List<ProtectEvidence> evidences = virusEvidence.evidence(testData);
        assertEquals(2, evidences.size());

        // The test data has reportable HPV virus
        assertTrue(findByEvent(evidences, "HPV Positive").reported());

        // The test data has non-reportable EBV virus
        assertFalse(findByEvent(evidences, "EBV Positive").reported());
    }

    @NotNull
    private static ProtectEvidence findByEvent(@NotNull List<ProtectEvidence> evidences, @NotNull String event) {
        for (ProtectEvidence evidence : evidences) {
            if (evidence.genomicEvent().equals(event)) {
                return evidence;
            }
        }

        throw new IllegalStateException("Could not find evidence with genomic event: " + event);
    }

    @NotNull
    private static VirusInterpreterData createTestVirusInterpreterData() {
        List<AnnotatedVirus> viruses = Lists.newArrayList();

        viruses.add(VirusTestFactory.testAnnotatedVirusBuilder().interpretation(VirusInterpretation.HPV).reported(true).build());
        viruses.add(VirusTestFactory.testAnnotatedVirusBuilder().interpretation(VirusInterpretation.EBV).reported(false).build());
        viruses.add(VirusTestFactory.testAnnotatedVirusBuilder().interpretation(VirusInterpretation.EBV).reported(false).build());
        viruses.add(VirusTestFactory.testAnnotatedVirusBuilder().interpretation(VirusInterpretation.MCV).reported(true).build());
        viruses.add(VirusTestFactory.testAnnotatedVirusBuilder().interpretation(null).reported(true).build());
        viruses.add(VirusTestFactory.testAnnotatedVirusBuilder().interpretation(null).reported(false).build());

        return ImmutableVirusInterpreterData.builder().allViruses(viruses).build();
    }

}