package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.purple.PurpleData;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristic;

import org.jetbrains.annotations.NotNull;

public class PurpleSignatureEvidence {

    @NotNull
    private final PersonalizedEvidenceFactory personalizedEvidenceFactory;
    @NotNull
    private final List<ActionableCharacteristic> actionableSignatures;

    public PurpleSignatureEvidence(@NotNull final PersonalizedEvidenceFactory personalizedEvidenceFactory,
            @NotNull final List<ActionableCharacteristic> actionableCharacteristics) {
        this.personalizedEvidenceFactory = personalizedEvidenceFactory;
        this.actionableSignatures = actionableCharacteristics.stream()
                .filter(x -> x.name() == TumorCharacteristic.MICROSATELLITE_UNSTABLE
                        || x.name() == TumorCharacteristic.HIGH_TUMOR_MUTATIONAL_LOAD)
                .collect(Collectors.toList());
    }

    @NotNull
    public List<ProtectEvidence> evidence(@NotNull PurpleData purpleData) {
        List<ProtectEvidence> result = Lists.newArrayList();
        for (ActionableCharacteristic signature : actionableSignatures) {
            switch (signature.name()) {
                case MICROSATELLITE_UNSTABLE: {
                    if (purpleData.microsatelliteStatus() == MicrosatelliteStatus.MSI) {
                        ProtectEvidence evidence = personalizedEvidenceFactory.somaticReportableEvidence(signature)
                                .genomicEvent("Microsatellite unstable")
                                .build();
                        result.add(evidence);
                    }
                    break;
                }
                case HIGH_TUMOR_MUTATIONAL_LOAD: {
                    if (purpleData.tumorMutationalLoadStatus() == TumorMutationalStatus.HIGH) {
                        ProtectEvidence evidence = personalizedEvidenceFactory.somaticReportableEvidence(signature)
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
