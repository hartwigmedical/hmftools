package com.hartwig.hmftools.ckb.interpretation.evidence;

import com.hartwig.hmftools.ckb.datamodelinterpretation.ImmutableCkbEntryInterpretation;
import com.hartwig.hmftools.ckb.datamodelinterpretation.indication.ImmutableIndication;
import com.hartwig.hmftools.ckb.interpretation.common.CommonInterpretationFactory;
import com.hartwig.hmftools.ckb.interpretation.common.therapyinterpretation.TherapyInterpretation;
import com.hartwig.hmftools.ckb.interpretation.common.therapyinterpretation.TherapyInterpretationFactory;
import com.hartwig.hmftools.ckb.interpretation.common.variantinterpretation.VariantInterpretationFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.IndicationInfo;
import com.hartwig.hmftools.ckb.json.common.TherapyInfo;
import com.hartwig.hmftools.ckb.json.indication.Indication;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.json.therapy.Therapy;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class EvidenceFactory {

    private static final Logger LOGGER = LogManager.getLogger(EvidenceFactory.class);

    private EvidenceFactory() {

    }

    public static void interpretEvidence(@NotNull MolecularProfile molecularProfile, @NotNull CkbJsonDatabase ckbEntry,
            @NotNull ImmutableCkbEntryInterpretation.Builder outputBuilder) {
        LOGGER.info(molecularProfile.id());
        for (EvidenceInfo evidenceInfo : molecularProfile.variantLevelEvidence().evidence()) {
            LOGGER.info(evidenceInfo.id());
            outputBuilder.addEvidenceInterpretation(ImmutableEvidenceInterpretation.builder()
                    .id(evidenceInfo.id())
                    .approvalStatus(evidenceInfo.approvalStatus())
                    .evidenceType(evidenceInfo.evidenceType())
                    .efficacyEvidence(evidenceInfo.efficacyEvidence())
                    .variantInterpretation(VariantInterpretationFactory.extractVariantGeneInfo(ckbEntry,
                            molecularProfile,
                            evidenceInfo.molecularProfile()).build())
                    .therapyInterpretation(extractTherapyEvidence(ckbEntry, evidenceInfo.therapy()))
                    .indications(extractIndication(ckbEntry, evidenceInfo.indication()))
                    .responseType(evidenceInfo.responseType())
                    .references(CommonInterpretationFactory.extractReferences(evidenceInfo.references(), ckbEntry))
                    .ampCapAscoEvidenceLevel(evidenceInfo.ampCapAscoEvidenceLevel())
                    .ampCapAscoInferredTier(evidenceInfo.ampCapAscoInferredTier())
                    .build());

        }
    }

    @Nullable
    private static TherapyInterpretation extractTherapyEvidence(@NotNull CkbJsonDatabase ckbEntry, @NotNull TherapyInfo therapyInfo) {
        TherapyInterpretation therapyInterpretation = null;
        for (Therapy therapy : ckbEntry.therapies()) {
            if (therapyInfo.id() == therapy.id()) {
                therapyInterpretation = TherapyInterpretationFactory.extractTherapyInterpretation(therapy, ckbEntry);
            }
        }
        return therapyInterpretation;
    }

    @NotNull
    private static com.hartwig.hmftools.ckb.datamodelinterpretation.indication.Indication extractIndication(
            @NotNull CkbJsonDatabase ckbEntry, @Nullable IndicationInfo indicationInfo) {
        ImmutableIndication.Builder outputBuilder = ImmutableIndication.builder();
        for (Indication indication : ckbEntry.indications()) {
            if (indicationInfo.id().equals(indication.id())) {
                outputBuilder
                        // TODO Switch to String for ID
                        .id(Integer.parseInt(indication.id()))
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
