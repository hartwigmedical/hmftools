package com.hartwig.hmftools.ckb.interpretation.evidence;

import com.hartwig.hmftools.ckb.datamodelinterpretation.indication.ImmutableIndication;
import com.hartwig.hmftools.ckb.interpretation.ImmutableCkbEntryInterpretation;
import com.hartwig.hmftools.ckb.interpretation.common.CommonInterpretationFactory;
import com.hartwig.hmftools.ckb.interpretation.common.therapyinterpretation.TherapyInterpretation;
import com.hartwig.hmftools.ckb.interpretation.common.therapyinterpretation.TherapyInterpretationFactory;
import com.hartwig.hmftools.ckb.interpretation.common.molecularprofileinterpretation.MolecularProfileInterpretationFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.IndicationInfo;
import com.hartwig.hmftools.ckb.json.common.TherapyInfo;
import com.hartwig.hmftools.ckb.json.indication.Indication;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.json.therapy.Therapy;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class EvidenceFactory {

    private EvidenceFactory() {

    }

    public static void interpretEvidence(@NotNull MolecularProfile molecularProfile, @NotNull CkbJsonDatabase ckbEntry,
            @NotNull ImmutableCkbEntryInterpretation.Builder outputBuilder) {
        for (EvidenceInfo evidenceInfo : molecularProfile.variantLevelEvidence().evidence()) {
            outputBuilder.addEvidenceInterpretations(ImmutableEvidenceInterpretation.builder()
                    .id(evidenceInfo.id())
                    .approvalStatus(evidenceInfo.approvalStatus())
                    .evidenceType(evidenceInfo.evidenceType())
                    .efficacyEvidence(evidenceInfo.efficacyEvidence())
                    .variantInterpretation(MolecularProfileInterpretationFactory.extractVariantGeneInfo(ckbEntry,
                            molecularProfile,
                            evidenceInfo.molecularProfile()).build())
                    .therapyInterpretation(extractTherapyEvidence(ckbEntry, evidenceInfo.therapy(), molecularProfile))
                    .indications(extractIndication(ckbEntry, evidenceInfo.indication()))
                    .responseType(evidenceInfo.responseType())
                    .references(CommonInterpretationFactory.extractReferences(evidenceInfo.references(), ckbEntry))
                    .ampCapAscoEvidenceLevel(evidenceInfo.ampCapAscoEvidenceLevel())
                    .ampCapAscoInferredTier(evidenceInfo.ampCapAscoInferredTier())
                    .build());

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

    @NotNull
    public static com.hartwig.hmftools.ckb.datamodelinterpretation.indication.Indication extractIndication(
            @NotNull CkbJsonDatabase ckbEntry, @Nullable IndicationInfo indicationInfo) {
        ImmutableIndication.Builder outputBuilder = ImmutableIndication.builder();
        for (Indication indication : ckbEntry.indications()) {
            if (indicationInfo.id().equals(indication.id())) {
                outputBuilder
                        .id(indication.id())
                        .name(indication.name())
                        .source(indication.source())
                        .definition(indication.definition())
                        .currentPreferredTerm(indication.currentPreferredTerm())
                        .lastUpdateDateFromDO(indication.lastUpdateDateFromDO())
                        .altIds(indication.altIds())
                        .termId(indication.termId());
            }
        }
        return outputBuilder.build();
    }
}