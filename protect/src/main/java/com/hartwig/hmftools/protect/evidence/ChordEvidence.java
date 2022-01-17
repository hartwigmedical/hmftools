package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristic;

import org.jetbrains.annotations.NotNull;

public class ChordEvidence {

    static final String HR_DEFICIENCY_EVENT = "HR deficiency";

    @NotNull
    private final PersonalizedEvidenceFactory personalizedEvidenceFactory;
    @NotNull
    private final List<ActionableCharacteristic> actionableCharacteristics;

    public ChordEvidence(@NotNull final PersonalizedEvidenceFactory personalizedEvidenceFactory,
            @NotNull final List<ActionableCharacteristic> actionableCharacteristics) {
        this.personalizedEvidenceFactory = personalizedEvidenceFactory;
        this.actionableCharacteristics = actionableCharacteristics.stream()
                .filter(x -> x.name() == TumorCharacteristic.HOMOLOGOUS_RECOMBINATION_DEFICIENT)
                .collect(Collectors.toList());
    }

    @NotNull
    public List<ProtectEvidence> evidence(@NotNull ChordAnalysis chordAnalysis) {
        List<ProtectEvidence> result = Lists.newArrayList();
        if (chordAnalysis.hrStatus() == ChordStatus.HR_DEFICIENT) {
            for (ActionableCharacteristic signature : actionableCharacteristics) {
                assert signature.name() == TumorCharacteristic.HOMOLOGOUS_RECOMBINATION_DEFICIENT;

                result.add(personalizedEvidenceFactory.somaticReportableEvidence(signature).event(HR_DEFICIENCY_EVENT).build());
            }
        }

        return result;
    }
}
