package com.hartwig.hmftools.ckb.datamodel.therapy;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.drug.DrugFactory;
import com.hartwig.hmftools.ckb.datamodel.reference.ReferenceFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.GlobalApprovalStatusInfo;
import com.hartwig.hmftools.ckb.json.common.TherapyInfo;
import com.hartwig.hmftools.ckb.json.therapy.JsonTherapy;

import org.jetbrains.annotations.NotNull;

public final class TherapyFactory {

    private TherapyFactory() {
    }

    @NotNull
    public static Therapy resolveTherapy(@NotNull CkbJsonDatabase ckbJsonDatabase, @NotNull TherapyInfo therapyInfo) {
        for (JsonTherapy therapy : ckbJsonDatabase.therapies()) {
            if (therapy.id() == therapyInfo.id()) {
                return ImmutableTherapy.builder()
                        .id(therapy.id())
                        .createDate(therapy.createDate())
                        .updateDate(therapy.updateDate())
                        .therapyName(therapy.therapyName())
                        .drugs(DrugFactory.extractDrugs(ckbJsonDatabase, therapy.drugs()))
                        .synonyms(therapy.synonyms())
                        .descriptions(extractTherapyDescriptions(ckbJsonDatabase, therapy.descriptions()))
                        .globalTherapyApprovalStatuses(convertGlobalApprovalStatuses(therapy.globalApprovalStatuses()))
                        .build();
            }
        }

        throw new IllegalStateException("Could not resolve CKB therapy with id '" + therapyInfo.id() + "'");
    }

    @NotNull
    private static List<TherapyDescription> extractTherapyDescriptions(@NotNull CkbJsonDatabase ckbJsonDatabase,
            @NotNull List<DescriptionInfo> descriptionInfos) {
        List<TherapyDescription> therapyDescriptions = Lists.newArrayList();

        for (DescriptionInfo descriptionInfo : descriptionInfos) {
            therapyDescriptions.add(ImmutableTherapyDescription.builder()
                    .description(descriptionInfo.description())
                    .references(ReferenceFactory.extractReferences(ckbJsonDatabase, descriptionInfo.references()))
                    .build());
        }
        return therapyDescriptions;
    }

    @NotNull
    private static List<GlobalTherapyApprovalStatus> convertGlobalApprovalStatuses(
            @NotNull List<GlobalApprovalStatusInfo> globalTherapyApprovalStatuses) {
        List<GlobalTherapyApprovalStatus> globalTherapyApprovalStatusesInterpretation = Lists.newArrayList();
        for (GlobalApprovalStatusInfo globalTherapyApprovalStatusInfo : globalTherapyApprovalStatuses) {
            globalTherapyApprovalStatusesInterpretation.add(ImmutableGlobalTherapyApprovalStatus.builder()
                    .id(globalTherapyApprovalStatusInfo.id())
                    .profileId(globalTherapyApprovalStatusInfo.molecularProfile().id())
                    .therapyId(globalTherapyApprovalStatusInfo.therapy().id())
                    .indicationId(globalTherapyApprovalStatusInfo.indication().id())
                    .approvalStatus(globalTherapyApprovalStatusInfo.approvalStatus())
                    .approvalAuthority(globalTherapyApprovalStatusInfo.approvalAuthority())
                    .build());
        }
        return globalTherapyApprovalStatusesInterpretation;
    }
}