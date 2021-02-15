package com.hartwig.hmftools.ckb.datamodel.therapy;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.drug.DrugFactory;
import com.hartwig.hmftools.ckb.datamodel.indication.IndicationFactory;
import com.hartwig.hmftools.ckb.datamodel.reference.ReferenceFactory;
import com.hartwig.hmftools.ckb.datamodel.variant.VariantFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.GlobalApprovalStatusInfo;
import com.hartwig.hmftools.ckb.json.molecularprofile.JsonMolecularProfile;
import com.hartwig.hmftools.ckb.json.therapy.JsonTherapy;

import org.jetbrains.annotations.NotNull;

public final class TherapyFactory {

    private TherapyFactory() {
    }

    @NotNull
    public static Therapy extractTherapy(@NotNull CkbJsonDatabase ckbJsonDatabase, @NotNull JsonTherapy therapy,
            @NotNull JsonMolecularProfile molecularProfile) {
        return ImmutableTherapy.builder()
                .id(therapy.id())
                .therapyName(therapy.therapyName())
                .synonyms(therapy.synonyms())
                .descriptions(extractTherapyDescriptions(ckbJsonDatabase, therapy.descriptions()))
                .createDate(therapy.createDate())
                .updateDate(therapy.updateDate())
                .drugs(DrugFactory.extractDrugs(ckbJsonDatabase, therapy.drugs()))
                .globalTherapyApprovalStatuses(extractGlobalApprovalStatuses(ckbJsonDatabase,
                        therapy.globalApprovalStatuses(),
                        molecularProfile,
                        therapy.id()))
                .build();
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
    private static List<GlobalTherapyApprovalStatus> extractGlobalApprovalStatuses(@NotNull CkbJsonDatabase ckbJsonDatabase,
            @NotNull List<GlobalApprovalStatusInfo> globalTherapyApprovalStatuses, @NotNull JsonMolecularProfile molecularProfile,
            int therapyId) {
        List<GlobalTherapyApprovalStatus> globalTherapyApprovalStatusesInterpretation = Lists.newArrayList();
        for (GlobalApprovalStatusInfo globalTherapyApprovalStatusInfo : globalTherapyApprovalStatuses) {
            if (therapyId == globalTherapyApprovalStatusInfo.therapy().id()
                    && molecularProfile.id() == globalTherapyApprovalStatusInfo.molecularProfile().id()) {
                globalTherapyApprovalStatusesInterpretation.add(ImmutableGlobalTherapyApprovalStatus.builder()
                        .id(globalTherapyApprovalStatusInfo.id())
                        .indication(IndicationFactory.extractIndication(ckbJsonDatabase, globalTherapyApprovalStatusInfo.indication()))
                        .variants(VariantFactory.extractVariants(ckbJsonDatabase, molecularProfile.geneVariants()))
                        .approvalStatus(globalTherapyApprovalStatusInfo.approvalStatus())
                        .approvalAuthority(globalTherapyApprovalStatusInfo.approvalAuthority())
                        .build());
            }
        }
        return globalTherapyApprovalStatusesInterpretation;
    }
}