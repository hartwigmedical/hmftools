package com.hartwig.hmftools.ckb.datamodel.evidence;

import com.hartwig.hmftools.ckb.datamodel.ImmutableCkbEntry;
import com.hartwig.hmftools.ckb.datamodel.common.CommonInterpretationFactory;
import com.hartwig.hmftools.ckb.datamodel.common.molecularprofile.MolecularProfileInterpretationFactory;
import com.hartwig.hmftools.ckb.datamodel.common.therapy.Therapy;
import com.hartwig.hmftools.ckb.datamodel.common.therapy.TherapyInterpretationFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.TherapyInfo;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfile;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class EvidenceFactory {

    private EvidenceFactory() {
    }

    public static void interpretVariantEvidence(@NotNull MolecularProfile molecularProfile, @NotNull CkbJsonDatabase ckbEntry,
            @NotNull ImmutableCkbEntry.Builder outputBuilder) {
        for (EvidenceInfo evidenceInfo : molecularProfile.variantLevelEvidence().evidence()) {
            if (molecularProfile.id() == evidenceInfo.molecularProfile().id()) {
                outputBuilder.addEvidences(ImmutableEvidence.builder()
                        .id(evidenceInfo.id())
                        .approvalStatus(evidenceInfo.approvalStatus())
                        .evidenceType(evidenceInfo.evidenceType())
                        .efficacyEvidence(evidenceInfo.efficacyEvidence())
                        .variants(MolecularProfileInterpretationFactory.extractVariantGeneInfo(ckbEntry,
                                molecularProfile))
                        .therapy(extractTherapyEvidence(ckbEntry, evidenceInfo.therapy(), molecularProfile))
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
    private static Therapy extractTherapyEvidence(@NotNull CkbJsonDatabase ckbEntry, @NotNull TherapyInfo therapyInfo,
            @NotNull MolecularProfile molecularProfile) {
        Therapy therapyInterpretation = null;
        for (com.hartwig.hmftools.ckb.json.therapy.Therapy therapy : ckbEntry.therapies()) {
            if (therapyInfo.id() == therapy.id()) {
                therapyInterpretation = TherapyInterpretationFactory.extractTherapy(therapy, ckbEntry, molecularProfile);
            }
        }
        return therapyInterpretation;
    }
}