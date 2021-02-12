package com.hartwig.hmftools.ckb.interpretation.common.therapyinterpretation;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodelinterpretation.globaltherapyapprovalstatus.GlobalTherapyApprovalStatus;
import com.hartwig.hmftools.ckb.datamodelinterpretation.globaltherapyapprovalstatus.ImmutableGlobalTherapyApprovalStatus;
import com.hartwig.hmftools.ckb.datamodelinterpretation.therapy.ImmutableTherapy;
import com.hartwig.hmftools.ckb.datamodelinterpretation.therapy.ImmutableTherapyDescription;
import com.hartwig.hmftools.ckb.datamodelinterpretation.therapy.TherapyDescription;
import com.hartwig.hmftools.ckb.interpretation.common.CommonInterpretationFactory;
import com.hartwig.hmftools.ckb.interpretation.common.druginterpretation.DrugInterpretationFactory;
import com.hartwig.hmftools.ckb.interpretation.common.molecularprofileinterpretation.MolecularProfileInterpretationFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.GlobalApprovalStatusInfo;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.json.therapy.Therapy;

import org.jetbrains.annotations.NotNull;

public class TherapyInterpretationFactory {

    private TherapyInterpretationFactory() {

    }

    @NotNull
    public static TherapyInterpretation extractTherapyInterpretation(@NotNull Therapy therapy, @NotNull CkbJsonDatabase ckbEntry,
            @NotNull MolecularProfile molecularProfile) {
        return ImmutableTherapyInterpretation.builder().therapy(extractTherapy(therapy, ckbEntry, molecularProfile)).build();
    }

    @NotNull
    private static com.hartwig.hmftools.ckb.datamodelinterpretation.therapy.Therapy extractTherapy(@NotNull Therapy therapy,
            @NotNull CkbJsonDatabase ckbEntry, @NotNull MolecularProfile molecularProfile) {
        return ImmutableTherapy.builder()
                .id(therapy.id())
                .therapyName(therapy.therapyName())
                .synonyms(therapy.synonyms())
                .descriptions(extractTherapyDescriptions(therapy.descriptions(), ckbEntry))
                .createDate(therapy.createDate())
                .updateDate(therapy.updateDate())
                .drugInterpretation(DrugInterpretationFactory.extractDrugInterpretation(therapy.drugs(), ckbEntry))
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
            if (therapyId == globalTherapyApprovalStatusInfo.therapy().id()) {
                globalTherapyApprovalStatusesInterpretation.add(ImmutableGlobalTherapyApprovalStatus.builder()
                        .id(globalTherapyApprovalStatusInfo.id())
                        .indications(CommonInterpretationFactory.extractIndication(ckbEntry, globalTherapyApprovalStatusInfo.indication()))
                        .variantInterpretation(MolecularProfileInterpretationFactory.extractVariantGeneInfo(ckbEntry,
                                molecularProfile,
                                globalTherapyApprovalStatusInfo.molecularProfile()).build())
                        .approvalStatus(globalTherapyApprovalStatusInfo.approvalStatus())
                        .approvalAuthority(globalTherapyApprovalStatusInfo.approvalAuthority())
                        .build());
            }
        }
        return globalTherapyApprovalStatusesInterpretation;
    }
}