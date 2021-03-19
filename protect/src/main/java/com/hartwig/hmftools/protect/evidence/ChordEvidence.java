package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableSignature;
import com.hartwig.hmftools.serve.extraction.characteristic.SignatureName;

import org.jetbrains.annotations.NotNull;

public class ChordEvidence {

    @NotNull
    private final PersonalizedEvidenceFactory personalizedEvidenceFactory;
    @NotNull
    private final List<ActionableSignature> actionableSignatures;

    public ChordEvidence(@NotNull final PersonalizedEvidenceFactory personalizedEvidenceFactory,
            @NotNull final List<ActionableSignature> actionableSignatures) {
        this.personalizedEvidenceFactory = personalizedEvidenceFactory;
        this.actionableSignatures = actionableSignatures.stream()
                .filter(x -> x.name() == SignatureName.HOMOLOGOUS_RECOMBINATION_DEFICIENT)
                .collect(Collectors.toList());
    }

    @NotNull
    public List<ProtectEvidence> evidence(@NotNull ChordAnalysis chordAnalysis) {
        List<ProtectEvidence> result = Lists.newArrayList();
        if (chordAnalysis.hrStatus() == ChordStatus.HR_DEFICIENT) {
            for (ActionableSignature signature : actionableSignatures) {
                assert signature.name() == SignatureName.HOMOLOGOUS_RECOMBINATION_DEFICIENT;

                result.add(personalizedEvidenceFactory.somaticallyReportableEvidence(signature)
                        .genomicEvent("HR deficiency signature")
                        .build());
            }
        }

        return result;
    }
}
