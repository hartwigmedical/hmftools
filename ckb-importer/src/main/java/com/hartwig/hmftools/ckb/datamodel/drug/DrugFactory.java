package com.hartwig.hmftools.ckb.datamodel.drug;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.reference.ReferenceFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.DrugClassInfo;
import com.hartwig.hmftools.ckb.json.common.DrugInfo;
import com.hartwig.hmftools.ckb.json.drug.JsonDrug;
import com.hartwig.hmftools.ckb.json.drugclass.JsonDrugClass;

import org.jetbrains.annotations.NotNull;

public final class DrugFactory {

    private DrugFactory() {
    }

    @NotNull
    public static List<Drug> extractDrugs(@NotNull CkbJsonDatabase ckbJsonDatabase, @NotNull List<DrugInfo> drugInfos) {
        List<Drug> drugs = Lists.newArrayList();
        for (DrugInfo drugInfo : drugInfos) {
            for (JsonDrug drug : ckbJsonDatabase.drugs()) {
                if (drugInfo.id() == drug.id()) {
                    drugs.add(ImmutableDrug.builder()
                            .id(drug.id())
                            .drugName(drug.drugName())
                            .terms(drug.terms())
                            .synonyms(drug.synonyms())
                            .tradeName(drug.tradeName())
                            .drugDescriptions(extractDrugDescriptions(ckbJsonDatabase, drug.descriptions()))
                            .casRegistryNum(drug.casRegistryNum())
                            .ncitId(drug.ncitId())
                            .createDate(drug.createDate())
                            .drugClasses(extractDrugClasses(ckbJsonDatabase, drug))
                            .build());
                }
            }
        }
        return drugs;
    }

    @NotNull
    private static List<DrugClass> extractDrugClasses(@NotNull CkbJsonDatabase ckbJsonDatabase, @NotNull JsonDrug drug) {
        List<DrugClass> drugClasses = Lists.newArrayList();
        for (DrugClassInfo drugClassInfo : drug.drugClasses()) {
            for (JsonDrugClass drugClass : ckbJsonDatabase.drugClasses()) {
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
    private static List<DrugDescription> extractDrugDescriptions(@NotNull CkbJsonDatabase ckbJsonDatabase,
            @NotNull List<DescriptionInfo> descriptionInfos) {
        List<DrugDescription> drugDescriptions = Lists.newArrayList();

        for (DescriptionInfo descriptionInfo : descriptionInfos) {
            drugDescriptions.add(ImmutableDrugDescription.builder()
                    .description(descriptionInfo.description())
                    .references(ReferenceFactory.extractReferences(ckbJsonDatabase, descriptionInfo.references()))
                    .build());
        }
        return drugDescriptions;
    }
}