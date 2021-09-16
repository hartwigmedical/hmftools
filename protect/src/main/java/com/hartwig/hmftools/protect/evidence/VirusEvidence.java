package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.virus.AnnotatedVirusV1;
import com.hartwig.hmftools.common.virus.VirusConstants;
import com.hartwig.hmftools.common.virus.VirusInterpreterData;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristic;

import org.jetbrains.annotations.NotNull;

public class VirusEvidence {

    @NotNull
    private final PersonalizedEvidenceFactory personalizedEvidenceFactory;
    @NotNull
    private final List<ActionableCharacteristic> actionableViruses;

    public VirusEvidence(@NotNull final PersonalizedEvidenceFactory personalizedEvidenceFactory,
            @NotNull final List<ActionableCharacteristic> actionableCharacteristics) {
        this.personalizedEvidenceFactory = personalizedEvidenceFactory;
        this.actionableViruses = actionableCharacteristics.stream()
                .filter(x -> x.name() == TumorCharacteristic.HPV_POSITIVE || x.name() == TumorCharacteristic.EBV_POSITIVE)
                .collect(Collectors.toList());
    }

    @NotNull
    public List<ProtectEvidence> evidence(@NotNull VirusInterpreterData virusInterpreterData) {
        List<AnnotatedVirusV1> hpv = virusesWithInterpretation(virusInterpreterData, VirusConstants.fromVirusName("HPV"));
        List<AnnotatedVirusV1> ebv = virusesWithInterpretation(virusInterpreterData, VirusConstants.fromVirusName("EBV"));

        List<ProtectEvidence> result = Lists.newArrayList();
        for (ActionableCharacteristic virus : actionableViruses) {
            switch (virus.name()) {
                case HPV_POSITIVE: {
                    if (!hpv.isEmpty()) {
                        ProtectEvidence evidence = personalizedEvidenceFactory.somaticEvidence(virus)
                                .genomicEvent("HPV Positive")
                                .reported(hpv.stream().anyMatch(x -> x.reported()))
                                .build();
                        result.add(evidence);
                    }
                    break;
                }
                case EBV_POSITIVE: {
                    if (!ebv.isEmpty()) {
                        ProtectEvidence evidence = personalizedEvidenceFactory.somaticEvidence(virus)
                                .genomicEvent("EBV Positive")
                                .reported(ebv.stream().anyMatch(x -> x.reported()))
                                .build();
                        result.add(evidence);
                    }
                    break;
                }
            }
        }
        return result;
    }

    @NotNull
    private static List<AnnotatedVirusV1> virusesWithInterpretation(@NotNull VirusInterpreterData virusInterpreterData,
            @NotNull VirusConstants interpretationToInclude) {
        List<AnnotatedVirusV1> virusesWithInterpretation = Lists.newArrayList();
        for (AnnotatedVirusV1 virus : virusInterpreterData.reportableViruses()) {
            if (virus.interpretation().equals(interpretationToInclude.name())) {
                virusesWithInterpretation.add(virus);
            }
        }

        for (AnnotatedVirusV1 virus : virusInterpreterData.unreportedViruses()) {
            if (virus.interpretation().equals(interpretationToInclude.name())) {
                virusesWithInterpretation.add(virus);
            }
        }
        return virusesWithInterpretation;
    }
}
