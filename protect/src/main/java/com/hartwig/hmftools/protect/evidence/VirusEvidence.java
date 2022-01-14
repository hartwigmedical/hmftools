package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceType;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusConstants;
import com.hartwig.hmftools.common.virus.VirusInterpreterData;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristic;

import org.jetbrains.annotations.NotNull;

public class VirusEvidence {

    static final String HPV_POSITIVE_EVENT = "HPV positive";
    static final String EBV_POSITIVE_EVENT = "EBV positive";

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
        List<AnnotatedVirus> hpv = virusesWithInterpretation(virusInterpreterData, VirusConstants.HPV);
        List<AnnotatedVirus> ebv = virusesWithInterpretation(virusInterpreterData, VirusConstants.EBV);

        List<ProtectEvidence> result = Lists.newArrayList();
        for (ActionableCharacteristic virus : actionableViruses) {
            switch (virus.name()) {
                case HPV_POSITIVE: {
                    if (!hpv.isEmpty()) {
                        ProtectEvidence evidence = personalizedEvidenceFactory.somaticEvidence(virus)
                                .reported(hpv.stream()
                                        .anyMatch(x -> x.reported() && x.virusDriverLikelihoodType().equals(VirusLikelihoodType.HIGH)))
                                .event(HPV_POSITIVE_EVENT)
                                .evidenceType(ProtectEvidenceType.VIRAL_PRESENCE)
                                .build();
                        result.add(evidence);
                    }
                    break;
                }
                case EBV_POSITIVE: {
                    if (!ebv.isEmpty()) {
                        ProtectEvidence evidence = personalizedEvidenceFactory.somaticEvidence(virus)
                                .reported(ebv.stream().anyMatch(x -> x.reported()))
                                .event(EBV_POSITIVE_EVENT)
                                .evidenceType(ProtectEvidenceType.VIRAL_PRESENCE)
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
    private static List<AnnotatedVirus> virusesWithInterpretation(@NotNull VirusInterpreterData virusInterpreterData,
            @NotNull VirusConstants interpretationToInclude) {
        List<AnnotatedVirus> virusesWithInterpretation = Lists.newArrayList();
        for (AnnotatedVirus virus : virusInterpreterData.reportableViruses()) {
            String interpretation = virus.interpretation();
            if (interpretation != null && interpretation.equals(interpretationToInclude.name())) {
                virusesWithInterpretation.add(virus);
            }
        }

        for (AnnotatedVirus virus : virusInterpreterData.unreportedViruses()) {
            String interpretation = virus.interpretation();
            if (interpretation != null && interpretation.equals(interpretationToInclude.name())) {
                virusesWithInterpretation.add(virus);
            }
        }
        return virusesWithInterpretation;
    }
}
