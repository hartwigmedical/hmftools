package com.hartwig.hmftools.ckb.datamodel.evidence;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.indication.IndicationFactory;
import com.hartwig.hmftools.ckb.datamodel.reference.ReferenceFactory;
import com.hartwig.hmftools.ckb.datamodel.therapy.TherapyFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class EvidenceFactory {

    private static final Logger LOGGER = LogManager.getLogger(EvidenceFactory.class);

    private EvidenceFactory() {
    }

    @NotNull
    public static List<Evidence> extractEvidences(@NotNull CkbJsonDatabase ckbJsonDatabase, @NotNull List<EvidenceInfo> evidenceInfos,
            int molecularProfileId) {
        List<Evidence> evidences = Lists.newArrayList();
        for (EvidenceInfo evidenceInfo : evidenceInfos) {
            if (evidenceInfo.molecularProfile().id() == molecularProfileId) {
                evidences.add(ImmutableEvidence.builder()
                        .id(evidenceInfo.id())
                        .therapy(TherapyFactory.resolveTherapy(ckbJsonDatabase, evidenceInfo.therapy()))
                        .indication(IndicationFactory.resolveIndication(ckbJsonDatabase, evidenceInfo.indication()))
                        .responseType(evidenceInfo.responseType())
                        .evidenceType(evidenceInfo.evidenceType())
                        .efficacyEvidence(evidenceInfo.efficacyEvidence())
                        .approvalStatus(evidenceInfo.approvalStatus())
                        .ampCapAscoEvidenceLevel(evidenceInfo.ampCapAscoEvidenceLevel())
                        .ampCapAscoInferredTier(evidenceInfo.ampCapAscoInferredTier())
                        .references(ReferenceFactory.extractReferences(ckbJsonDatabase, evidenceInfo.references()))
                        .build());
            } else {
                LOGGER.warn("Variant level evidence on profile {} for different profile: {}",
                        molecularProfileId,
                        evidenceInfo.molecularProfile().id());
            }
        }
        return evidences;
    }
}