package com.hartwig.hmftools.ckb.datamodel.common.drug;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.common.CommonInterpretationFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.DrugClassInfo;
import com.hartwig.hmftools.ckb.json.common.DrugInfo;

import org.jetbrains.annotations.NotNull;

public final class DrugInterpretationFactory {

    private DrugInterpretationFactory() {
    }

    @NotNull
    public static List<Drug> extractDrugs(@NotNull List<DrugInfo> drugInfos, @NotNull CkbJsonDatabase ckbEntry) {
        List<Drug> drugs = Lists.newArrayList();
        for (DrugInfo drugInfo : drugInfos) {
            for (com.hartwig.hmftools.ckb.json.drug.Drug drug : ckbEntry.drugs()) {
                if (drugInfo.id() == drug.id()) {
                    drugs.add(ImmutableDrug.builder()
                            .id(drug.id())
                            .drugName(drug.drugName())
                            .terms(drug.terms())
                            .synonyms(drug.synonyms())
                            .tradeName(drug.tradeName())
                            .drugDescriptions(extractDrugDescriptions(drug.descriptions(), ckbEntry))
                            .casRegistryNum(drug.casRegistryNum())
                            .ncitId(drug.ncitId())
                            .createDate(drug.createDate())
                            .drugClasses(extractDrugClasses(drug, ckbEntry))
                            .build());
                }
            }
        }
        return drugs;
    }

    @NotNull
    private static List<DrugClass> extractDrugClasses(@NotNull com.hartwig.hmftools.ckb.json.drug.Drug drug,
            @NotNull CkbJsonDatabase ckbEntry) {
        List<com.hartwig.hmftools.ckb.datamodel.common.drug.DrugClass> drugClasses = Lists.newArrayList();
        for (DrugClassInfo drugClassInfo : drug.drugClasses()) {
            for (com.hartwig.hmftools.ckb.json.drugclass.DrugClass drugClass : ckbEntry.drugClasses()) {
                if (drugClassInfo.id() == drugClass.id()) {
                    drugClasses.add(ImmutableDrugClass.builder()
                            .id(drugClass.id())
                            .drugClass(drugClass.drugClass())
                            .createDate(drugClass.createDate())
                            .build());
                }
            }
        }
        return drugClasses;
    }

    @NotNull
    private static List<DrugDescription> extractDrugDescriptions(@NotNull List<DescriptionInfo> descriptionInfos,
            @NotNull CkbJsonDatabase ckbEntry) {
        List<DrugDescription> drugDescriptions = Lists.newArrayList();

        for (DescriptionInfo descriptionInfo : descriptionInfos) {
            drugDescriptions.add(ImmutableDrugDescription.builder()
                    .description(descriptionInfo.description())
                    .references(CommonInterpretationFactory.extractReferences(descriptionInfo.references(), ckbEntry))
                    .build());
        }
        return drugDescriptions;
    }
}