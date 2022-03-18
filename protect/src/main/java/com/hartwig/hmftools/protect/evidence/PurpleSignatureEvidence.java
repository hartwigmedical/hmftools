package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.purple.PurpleData;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicAnnotation;

import org.jetbrains.annotations.NotNull;

public class PurpleSignatureEvidence {

    static final String MICROSATELLITE_UNSTABLE_EVENT = "Microsatellite unstable";
    static final String HIGH_TUMOR_LOAD_EVENT = "High tumor mutation load";
    static final String HIGH_TUMOR_BURDEN_EVENT = "High tumor mutation burden";

    @NotNull
    private final PersonalizedEvidenceFactory personalizedEvidenceFactory;
    @NotNull
    private final List<ActionableCharacteristic> actionableSignatures;

    public PurpleSignatureEvidence(@NotNull final PersonalizedEvidenceFactory personalizedEvidenceFactory,
            @NotNull final List<ActionableCharacteristic> actionableCharacteristics) {
        this.personalizedEvidenceFactory = personalizedEvidenceFactory;
        this.actionableSignatures = actionableCharacteristics.stream()
                .filter(x -> x.name() == TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE
                        || x.name() == TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD
                        || x.name() == TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_BURDEN)
                .collect(Collectors.toList());
    }

    @NotNull
    public List<ProtectEvidence> evidence(@NotNull PurpleData purpleData) {
        List<ProtectEvidence> result = Lists.newArrayList();
        for (ActionableCharacteristic signatureMSI : actionableSignatures) {
            if (signatureMSI.name() == TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE) {
                if (purpleData.microsatelliteStatus() == MicrosatelliteStatus.MSI) {
                    switch (signatureMSI.comparator()) {
                        case GREATER:
                            if (purpleData.microsatelliteIndelsPerMb() > signatureMSI.cutoff()) {
                                result.add(generateMSIEvidences(signatureMSI));
                            }
                            break;
                        case EQUAL_OR_GREATER:
                            if (purpleData.microsatelliteIndelsPerMb() >= signatureMSI.cutoff()) {
                                result.add(generateMSIEvidences(signatureMSI));
                            }
                            break;
                        case LOWER:
                            if (purpleData.microsatelliteIndelsPerMb() < signatureMSI.cutoff()) {
                                result.add(generateMSIEvidences(signatureMSI));
                            }
                            break;
                        case EQUAL_OR_LOWER:
                            if (purpleData.microsatelliteIndelsPerMb() <= signatureMSI.cutoff()) {
                                result.add(generateMSIEvidences(signatureMSI));
                            }
                            break;
                    }
                }
            }
        }

        for (ActionableCharacteristic signatureTML : actionableSignatures) {
            if (signatureTML.name() == TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD) {
                if (purpleData.tumorMutationalLoadStatus() == TumorMutationalStatus.HIGH) {
                    switch (signatureTML.comparator()) {
                        case EQUAL_OR_LOWER:
                            if (purpleData.tumorMutationalLoad() <= signatureTML.cutoff()) {
                                result.add(generateMTLEvidences(signatureTML));
                            }
                            break;
                        case LOWER:
                            if (purpleData.tumorMutationalLoad() < signatureTML.cutoff()) {
                                result.add(generateMTLEvidences(signatureTML));
                            }
                            break;
                        case EQUAL_OR_GREATER:
                            if (purpleData.tumorMutationalLoad() >= signatureTML.cutoff()) {
                                result.add(generateMTLEvidences(signatureTML));
                            }
                            break;
                        case GREATER:
                            if (purpleData.tumorMutationalLoad() > signatureTML.cutoff()) {
                                result.add(generateMTLEvidences(signatureTML));
                            }
                            break;
                    }
                }
            }
        }

        for (ActionableCharacteristic signatureTMB : actionableSignatures) {
            if (signatureTMB.name() == TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_BURDEN) {
                switch (signatureTMB.comparator()) {
                    case EQUAL_OR_LOWER:
                        if (purpleData.tumorMutationalBurdenPerMb() <= signatureTMB.cutoff()) {
                            result.add(generateMTBEvidences(signatureTMB));
                        }
                        break;
                    case LOWER:
                        if (purpleData.tumorMutationalBurdenPerMb() < signatureTMB.cutoff()) {
                            result.add(generateMTBEvidences(signatureTMB));
                        }
                        break;
                    case EQUAL_OR_GREATER:
                        if (purpleData.tumorMutationalBurdenPerMb() >= signatureTMB.cutoff()) {
                            result.add(generateMTBEvidences(signatureTMB));
                        }
                        break;
                    case GREATER:
                        if (purpleData.tumorMutationalBurdenPerMb() > signatureTMB.cutoff()) {
                            result.add(generateMTBEvidences(signatureTMB));
                        }
                        break;
                }
            }
        }
        return result;
    }

    @NotNull
    public ProtectEvidence generateMSIEvidences(@NotNull ActionableCharacteristic signature) {
        return personalizedEvidenceFactory.somaticReportableEvidence(signature)
                .event(MICROSATELLITE_UNSTABLE_EVENT)
                .eventIsHighDriver(EvidenceDriverLikelihood.interpretSignatures())
                .build();
    }

    @NotNull
    public ProtectEvidence generateMTLEvidences(@NotNull ActionableCharacteristic signature) {
        return personalizedEvidenceFactory.somaticReportableEvidence(signature)
                .event(HIGH_TUMOR_LOAD_EVENT)
                .eventIsHighDriver(EvidenceDriverLikelihood.interpretSignatures())
                .build();
    }

    @NotNull
    public ProtectEvidence generateMTBEvidences(@NotNull ActionableCharacteristic signature) {
        return personalizedEvidenceFactory.somaticReportableEvidence(signature)
                .event(HIGH_TUMOR_BURDEN_EVENT)
                .eventIsHighDriver(EvidenceDriverLikelihood.interpretSignatures())
                .build();
    }
}
