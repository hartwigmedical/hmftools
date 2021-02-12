package com.hartwig.hmftools.ckb.interpretation.evidence;

import com.hartwig.hmftools.ckb.interpretation.ImmutableCkbEntryInterpretation;
import com.hartwig.hmftools.ckb.interpretation.common.CommonInterpretationFactory;
import com.hartwig.hmftools.ckb.interpretation.common.therapyinterpretation.TherapyInterpretation;
import com.hartwig.hmftools.ckb.interpretation.common.therapyinterpretation.TherapyInterpretationFactory;
import com.hartwig.hmftools.ckb.interpretation.common.molecularprofileinterpretation.MolecularProfileInterpretationFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.TherapyInfo;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.json.therapy.Therapy;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class EvidenceFactory {

    private EvidenceFactory() {

    }

    public static void interpretVariantEvidence(@NotNull MolecularProfile molecularProfile, @NotNull CkbJsonDatabase ckbEntry,
            @NotNull ImmutableCkbEntryInterpretation.Builder outputBuilder) {
        for (EvidenceInfo evidenceInfo : molecularProfile.variantLevelEvidence().evidence()) {
            if (molecularProfile.id() == evidenceInfo.molecularProfile().id()) {
                outputBuilder.addEvidenceInterpretations(ImmutableEvidenceInterpretation.builder()
                        .id(evidenceInfo.id())
                        .approvalStatus(evidenceInfo.approvalStatus())
                        .evidenceType(evidenceInfo.evidenceType())
                        .efficacyEvidence(evidenceInfo.efficacyEvidence())
                        .variantInterpretation(MolecularProfileInterpretationFactory.extractVariantGeneInfo(ckbEntry,
                                molecularProfile).build())
                        .therapyInterpretation(extractTherapyEvidence(ckbEntry, evidenceInfo.therapy(), molecularProfile))
                        .indication(CommonInterpretationFactory.extractIndication(ckbEntry, evidenceInfo.indication()))
                        .responseType(evidenceInfo.responseType())
                        .references(CommonInterpretationFactory.extractReferences(evidenceInfo.references(), ckbEntry))
                        .ampCapAscoEvidenceLevel(evidenceInfo.ampCapAscoEvidenceLevel())
                        .ampCapAscoInferredTier(evidenceInfo.ampCapAscoInferredTier())
                        .build());
            }
        }
    }

    @Nullable
    private static TherapyInterpretation extractTherapyEvidence(@NotNull CkbJsonDatabase ckbEntry, @NotNull TherapyInfo therapyInfo,
            @NotNull MolecularProfile molecularProfile) {
        TherapyInterpretation therapyInterpretation = null;
        for (Therapy therapy : ckbEntry.therapies()) {
            if (therapyInfo.id() == therapy.id()) {
                therapyInterpretation = TherapyInterpretationFactory.extractTherapyInterpretation(therapy, ckbEntry, molecularProfile);
            }
        }
        return therapyInterpretation;
    }
}