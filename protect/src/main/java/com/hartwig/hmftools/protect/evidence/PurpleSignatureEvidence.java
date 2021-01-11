package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.protect.purple.PurpleData;
import com.hartwig.hmftools.serve.actionability.signature.ActionableSignature;
import com.hartwig.hmftools.serve.extraction.signature.SignatureName;

import org.jetbrains.annotations.NotNull;

public class PurpleSignatureEvidence {

    @NotNull
    private final List<ActionableSignature> actionableSignatures;

    public PurpleSignatureEvidence(@NotNull final List<ActionableSignature> actionableSignatures) {
        this.actionableSignatures = actionableSignatures.stream()
                .filter(x -> x.name() == SignatureName.MICROSATELLITE_UNSTABLE || x.name() == SignatureName.HIGH_TUMOR_MUTATIONAL_LOAD)
                .collect(Collectors.toList());
    }

    @NotNull
    public List<ProtectEvidence> evidence(@NotNull Set<String> doids, @NotNull PurpleData purpleData) {
        List<ProtectEvidence> result = Lists.newArrayList();
        for (ActionableSignature signature : actionableSignatures) {
            ImmutableProtectEvidence.Builder signatureEvidenceBuilder =
                    ProtectEvidenceFunctions.builder(doids, signature).germline(false).reported(true);
            switch (signature.name()) {
                case MICROSATELLITE_UNSTABLE: {
                    if (purpleData.microsatelliteStatus() == MicrosatelliteStatus.MSI) {
                        ProtectEvidence evidence = signatureEvidenceBuilder.genomicEvent(MicrosatelliteStatus.MSI.display()).build();
                        result.add(evidence);
                    }
                    break;
                }
                case HIGH_TUMOR_MUTATIONAL_LOAD: {
                    if (purpleData.tumorMutationalLoadStatus() == TumorMutationalStatus.HIGH) {
                        ProtectEvidence evidence = signatureEvidenceBuilder.genomicEvent(TumorMutationalStatus.HIGH.display()).build();
                        result.add(evidence);
                    }
                    break;
                }
            }
        }
        return result;
    }
}
