package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.protect.characteristic.CharacteristicsFunctions;
import com.hartwig.serve.datamodel.characteristic.ActionableCharacteristic;
import com.hartwig.serve.datamodel.characteristic.TumorCharacteristicAnnotation;

import org.jetbrains.annotations.NotNull;

public class ChordEvidence {

    static final String HR_DEFICIENCY_EVENT = "HR deficiency";

    @NotNull
    private final PersonalizedEvidenceFactory personalizedEvidenceFactory;
    @NotNull
    private final List<ActionableCharacteristic> actionableCharacteristics;

    public ChordEvidence(@NotNull final PersonalizedEvidenceFactory personalizedEvidenceFactory,
            @NotNull final List<ActionableCharacteristic> actionableCharacteristics) {
        this.personalizedEvidenceFactory = personalizedEvidenceFactory;
        this.actionableCharacteristics = actionableCharacteristics.stream()
                .filter(x -> x.name() == TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT)
                .collect(Collectors.toList());
    }

    @NotNull
    public List<ProtectEvidence> evidence(@NotNull ChordData chordAnalysis) {
        List<ProtectEvidence> result = Lists.newArrayList();
        for (ActionableCharacteristic characteristic : actionableCharacteristics) {
            if (CharacteristicsFunctions.hasExplicitCutoff(characteristic)) {
                if (CharacteristicsFunctions.evaluateVersusCutoff(characteristic, chordAnalysis.hrdValue())) {
                    result.add(toHRDEvidence(characteristic));
                }
            } else if (chordAnalysis.hrStatus() == ChordStatus.HR_DEFICIENT) {
                result.add(toHRDEvidence(characteristic));
            }
        }

        return result;
    }

    @NotNull
    private ProtectEvidence toHRDEvidence(@NotNull ActionableCharacteristic signature) {
        return personalizedEvidenceFactory.somaticReportableEvidence(signature).event(HR_DEFICIENCY_EVENT).eventIsHighDriver(null).build();
    }
}
