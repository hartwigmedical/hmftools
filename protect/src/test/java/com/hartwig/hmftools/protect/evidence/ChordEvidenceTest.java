package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.chord.ChordTestFactory;
import com.hartwig.hmftools.common.chord.ImmutableChordAnalysis;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.serve.ServeTestFactory;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.actionability.characteristic.ImmutableActionableCharacteristic;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicAnnotation;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ChordEvidenceTest {

    @Test
    public void canDetermineEvidenceForChord() {
        ActionableCharacteristic signature1 = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic())
                .name(TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT)
                .build();

        ActionableCharacteristic signature2 = ImmutableActionableCharacteristic.builder()
                .from(ServeTestFactory.createTestActionableCharacteristic())
                .name(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD)
                .build();

        ChordEvidence chordEvidence =
                new ChordEvidence(EvidenceTestFactory.createTestEvidenceFactory(), Lists.newArrayList(signature1, signature2));

        ChordAnalysis hrDeficient = chordAnalysisWithStatus(ChordStatus.HR_DEFICIENT);
        List<ProtectEvidence> evidence = chordEvidence.evidence(hrDeficient);

        assertEquals(1, evidence.size());
        assertTrue(evidence.get(0).reported());
        assertEquals(ChordEvidence.HR_DEFICIENCY_EVENT, evidence.get(0).event());

        ChordAnalysis hrProficient = chordAnalysisWithStatus(ChordStatus.HR_PROFICIENT);
        assertTrue(chordEvidence.evidence(hrProficient).isEmpty());
    }

    @NotNull
    private static ChordAnalysis chordAnalysisWithStatus(@NotNull ChordStatus status) {
        return ImmutableChordAnalysis.builder().from(ChordTestFactory.createMinimalTestChordAnalysis()).hrStatus(status).build();
    }
}