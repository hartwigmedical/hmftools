package com.hartwig.hmftools.protect.evidence;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.lilac.LilacAllele;
import com.hartwig.hmftools.common.lilac.LilacData;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.serve.actionability.immuno.ActionableHLA;

import org.jetbrains.annotations.NotNull;

public class HlaEvidence {

    @NotNull
    private final PersonalizedEvidenceFactory personalizedEvidenceFactory;
    @NotNull
    private final List<ActionableHLA> actionableHLA;

    public HlaEvidence(@NotNull final PersonalizedEvidenceFactory personalizedEvidenceFactory,
            @NotNull final List<ActionableHLA> actionableHLA) {
        this.personalizedEvidenceFactory = personalizedEvidenceFactory;
        this.actionableHLA = actionableHLA;
    }

    public List<ProtectEvidence> evidence(@NotNull LilacData lilacData) {
        List<ProtectEvidence> result = Lists.newArrayList();
        for (LilacAllele lilacAllele : lilacData.alleles()) {
            result.addAll(evidence(lilacAllele, lilacData.qc()));
        }

        return result;
    }

    @NotNull
    private List<ProtectEvidence> evidence(@NotNull LilacAllele lilacAllele, @NotNull String lilacQc) {
        List<ProtectEvidence> result = Lists.newArrayList();

        for (ActionableHLA hla : actionableHLA) {
            if (hla.hlaType().equals(lilacAllele.name())) {
                    ProtectEvidence evidence = personalizedEvidenceFactory.evidenceBuilder(hla)
                            .gene(hla.hlaType())
                            .isCanonical(false)
                            .event(lilacAllele.name())
                            .reported(lilacQc.equals("PASS"))
                            .germline(false)
                            .build();
                    result.add(evidence);
            }
        }
        return result;
    }
}