package com.hartwig.hmftools.ckb.datamodel.common.therapy;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.common.CommonInterpretationFactory;
import com.hartwig.hmftools.ckb.datamodel.common.drug.DrugInterpretationFactory;
import com.hartwig.hmftools.ckb.datamodel.common.molecularprofile.MolecularProfileInterpretationFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.GlobalApprovalStatusInfo;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfile;

import org.jetbrains.annotations.NotNull;

public final class TherapyInterpretationFactory {

    private TherapyInterpretationFactory() {
    }

    @NotNull
    public static Therapy extractTherapy(@NotNull com.hartwig.hmftools.ckb.json.therapy.Therapy therapy, @NotNull CkbJsonDatabase ckbEntry,
            @NotNull MolecularProfile molecularProfile) {
        return ImmutableTherapy.builder()
                .id(therapy.id())
                .therapyName(therapy.therapyName())
                .synonyms(therapy.synonyms())
                .descriptions(extractTherapyDescriptions(therapy.descriptions(), ckbEntry))
                .createDate(therapy.createDate())
                .updateDate(therapy.updateDate())
                .drugs(DrugInterpretationFactory.extractDrugs(therapy.drugs(), ckbEntry))
                .globalTherapyApprovalStatuses(extractGlobalApprovalStatuses(therapy.globalApprovalStatuses(),
                        ckbEntry,
                        molecularProfile,
                        therapy.id()))
                .build();
    }

    @NotNull
    private static List<TherapyDescription> extractTherapyDescriptions(@NotNull List<DescriptionInfo> descriptionInfos,
            @NotNull CkbJsonDatabase ckbEntry) {
        List<TherapyDescription> therapyDescriptions = Lists.newArrayList();

        for (DescriptionInfo descriptionInfo : descriptionInfos) {
            therapyDescriptions.add(ImmutableTherapyDescription.builder()
                    .description(descriptionInfo.description())
                    .references(CommonInterpretationFactory.extractReferences(descriptionInfo.references(), ckbEntry))
                    .build());
        }
        return therapyDescriptions;
    }

    @NotNull
    private static List<GlobalTherapyApprovalStatus> extractGlobalApprovalStatuses(
            @NotNull List<GlobalApprovalStatusInfo> globalTherapyApprovalStatuses, @NotNull CkbJsonDatabase ckbEntry,
            @NotNull MolecularProfile molecularProfile, int therapyId) {
        List<GlobalTherapyApprovalStatus> globalTherapyApprovalStatusesInterpretation = Lists.newArrayList();
        for (GlobalApprovalStatusInfo globalTherapyApprovalStatusInfo : globalTherapyApprovalStatuses) {
            if (therapyId == globalTherapyApprovalStatusInfo.therapy().id()
                    && molecularProfile.id() == globalTherapyApprovalStatusInfo.molecularProfile().id()) {
                globalTherapyApprovalStatusesInterpretation.add(ImmutableGlobalTherapyApprovalStatus.builder()
                        .id(globalTherapyApprovalStatusInfo.id())
                        .indication(CommonInterpretationFactory.extractIndication(ckbEntry, globalTherapyApprovalStatusInfo.indication()))
                        .variants(MolecularProfileInterpretationFactory.extractVariantGeneInfo(ckbEntry, molecularProfile))
                        .approvalStatus(globalTherapyApprovalStatusInfo.approvalStatus())
                        .approvalAuthority(globalTherapyApprovalStatusInfo.approvalAuthority())
                        .build());
            }
        }
        return globalTherapyApprovalStatusesInterpretation;
    }
}