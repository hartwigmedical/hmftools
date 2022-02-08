package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusConstants;
import com.hartwig.hmftools.common.virus.VirusInterpreterData;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicAnnotation;

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
                .filter(x -> x.name() == TumorCharacteristicAnnotation.HPV_POSITIVE
                        || x.name() == TumorCharacteristicAnnotation.EBV_POSITIVE)
                .collect(Collectors.toList());
    }

    @NotNull
    public List<ProtectEvidence> evidence(@NotNull VirusInterpreterData virusInterpreterData) {
        List<AnnotatedVirus> hpv = virusesWithInterpretation(virusInterpreterData, VirusConstants.HPV);
        List<AnnotatedVirus> ebv = virusesWithInterpretation(virusInterpreterData, VirusConstants.EBV);

        boolean reportHPV = hasReportedWithHighDriverLikelihood(hpv);
        boolean reportEBV = hasReportedWithHighDriverLikelihood(ebv);

        List<ProtectEvidence> result = Lists.newArrayList();
        for (ActionableCharacteristic virus : actionableViruses) {
            switch (virus.name()) {
                case HPV_POSITIVE: {
                    if (!hpv.isEmpty()) {
                        ProtectEvidence evidence =
                                personalizedEvidenceFactory.somaticEvidence(virus).reported(reportHPV).event(HPV_POSITIVE_EVENT).build();
                        result.add(evidence);
                    }
                    break;
                }
                case EBV_POSITIVE: {
                    if (!ebv.isEmpty()) {
                        ProtectEvidence evidence =
                                personalizedEvidenceFactory.somaticEvidence(virus).reported(reportEBV).event(EBV_POSITIVE_EVENT).build();
                        result.add(evidence);
                    }
                    break;
                }
            }
        }
        return result;
    }

    private static boolean hasReportedWithHighDriverLikelihood(@NotNull List<AnnotatedVirus> viruses) {
        for (AnnotatedVirus virus : viruses) {
            if (virus.reported() && virus.virusDriverLikelihoodType() == VirusLikelihoodType.HIGH) {
                return true;
            }
        }

        return false;
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
