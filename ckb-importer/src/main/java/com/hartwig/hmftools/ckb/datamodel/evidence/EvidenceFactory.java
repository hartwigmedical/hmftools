package com.hartwig.hmftools.ckb.datamodel.evidence;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.common.CommonInterpretationFactory;
import com.hartwig.hmftools.ckb.datamodel.common.therapy.Therapy;
import com.hartwig.hmftools.ckb.datamodel.common.therapy.TherapyFactory;
import com.hartwig.hmftools.ckb.datamodel.common.variant.VariantFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.TherapyInfo;
import com.hartwig.hmftools.ckb.json.molecularprofile.JsonMolecularProfile;
import com.hartwig.hmftools.ckb.json.therapy.JsonTherapy;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class EvidenceFactory {

    private EvidenceFactory() {
    }

    @NotNull
    public static List<Evidence> interpretVariantEvidence(@NotNull CkbJsonDatabase ckbJsonDatabase,
            @NotNull JsonMolecularProfile molecularProfile) {
        List<Evidence> evidences = Lists.newArrayList();
        for (EvidenceInfo evidenceInfo : molecularProfile.variantLevelEvidence().evidences()) {
            if (molecularProfile.id() == evidenceInfo.molecularProfile().id()) {
                evidences.add(ImmutableEvidence.builder()
                        .id(evidenceInfo.id())
                        .approvalStatus(evidenceInfo.approvalStatus())
                        .evidenceType(evidenceInfo.evidenceType())
                        .efficacyEvidence(evidenceInfo.efficacyEvidence())
                        .variants(VariantFactory.extractVariants(ckbJsonDatabase, molecularProfile.geneVariants()))
                        .therapy(extractTherapyEvidence(ckbJsonDatabase, evidenceInfo.therapy(), molecularProfile))
                        .indication(CommonInterpretationFactory.extractIndication(ckbJsonDatabase, evidenceInfo.indication()))
                        .responseType(evidenceInfo.responseType())
                        .references(CommonInterpretationFactory.extractReferences(ckbJsonDatabase, evidenceInfo.references()))
                        .ampCapAscoEvidenceLevel(evidenceInfo.ampCapAscoEvidenceLevel())
                        .ampCapAscoInferredTier(evidenceInfo.ampCapAscoInferredTier())
                        .build());
            }
        }
        return evidences;
    }

    @Nullable
    private static Therapy extractTherapyEvidence(@NotNull CkbJsonDatabase ckbJsonDatabase, @NotNull TherapyInfo therapyInfo,
            @NotNull JsonMolecularProfile molecularProfile) {
        Therapy therapyInterpretation = null;
        for (JsonTherapy therapy : ckbJsonDatabase.therapies()) {
            if (therapyInfo.id() == therapy.id()) {
                therapyInterpretation = TherapyFactory.extractTherapy(ckbJsonDatabase, therapy, molecularProfile);
            }
        }
        return therapyInterpretation;
    }
}