package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.serve.actionability.signature.ActionableSignature;
import com.hartwig.hmftools.serve.extraction.signature.SignatureName;

import org.jetbrains.annotations.NotNull;

public class ChordEvidence {

    @NotNull
    private final List<ActionableSignature> actionableSignatures;

    public ChordEvidence(@NotNull final List<ActionableSignature> actionableSignatures) {
        this.actionableSignatures = actionableSignatures.stream()
                .filter(x -> x.name() == SignatureName.HOMOLOGOUS_RECOMBINATION_DEFICIENT)
                .collect(Collectors.toList());
    }

    @NotNull
    public List<ProtectEvidence> evidence(@NotNull Set<String> doids, @NotNull ChordAnalysis chordAnalysis) {
        List<ProtectEvidence> result = Lists.newArrayList();
        if (chordAnalysis.hrStatus() == ChordStatus.HR_DEFICIENT) {
            for (ActionableSignature signature : actionableSignatures) {
                assert signature.name() == SignatureName.HOMOLOGOUS_RECOMBINATION_DEFICIENT;

                result.add(ProtectEvidenceFunctions.builder(doids, signature)
                        .genomicEvent(ChordStatus.HR_DEFICIENT.display())
                        .germline(false)
                        .reported(true)
                        .build());
            }
        }
        return result;
    }
}
