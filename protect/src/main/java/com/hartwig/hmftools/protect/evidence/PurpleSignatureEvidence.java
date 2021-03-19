package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.protect.purple.PurpleData;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableSignature;
import com.hartwig.hmftools.serve.extraction.signature.SignatureName;

import org.jetbrains.annotations.NotNull;

public class PurpleSignatureEvidence {

    @NotNull
    private final PersonalizedEvidenceFactory personalizedEvidenceFactory;
    @NotNull
    private final List<ActionableSignature> actionableSignatures;

    public PurpleSignatureEvidence(@NotNull final PersonalizedEvidenceFactory personalizedEvidenceFactory,
            @NotNull final List<ActionableSignature> actionableSignatures) {
        this.personalizedEvidenceFactory = personalizedEvidenceFactory;
        this.actionableSignatures = actionableSignatures.stream()
                .filter(x -> x.name() == SignatureName.MICROSATELLITE_UNSTABLE || x.name() == SignatureName.HIGH_TUMOR_MUTATIONAL_LOAD)
                .collect(Collectors.toList());
    }

    @NotNull
    public List<ProtectEvidence> evidence(@NotNull PurpleData purpleData) {
        List<ProtectEvidence> result = Lists.newArrayList();
        for (ActionableSignature signature : actionableSignatures) {
            switch (signature.name()) {
                case MICROSATELLITE_UNSTABLE: {
                    if (purpleData.microsatelliteStatus() == MicrosatelliteStatus.MSI) {
                        ProtectEvidence evidence = personalizedEvidenceFactory.somaticallyReportableEvidence(signature)
                                .genomicEvent("Microsatellite unstable")
                                .build();
                        result.add(evidence);
                    }
                    break;
                }
                case HIGH_TUMOR_MUTATIONAL_LOAD: {
                    if (purpleData.tumorMutationalLoadStatus() == TumorMutationalStatus.HIGH) {
                        ProtectEvidence evidence = personalizedEvidenceFactory.somaticallyReportableEvidence(signature)
                                .genomicEvent("High tumor mutation load")
                                .build();
                        result.add(evidence);
                    }
                    break;
                }
            }
        }

        return result;
    }
}
