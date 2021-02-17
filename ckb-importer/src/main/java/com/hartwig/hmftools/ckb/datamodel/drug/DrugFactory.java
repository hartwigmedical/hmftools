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
            drugs.add(resolveDrug(ckbJsonDatabase, drugInfo));
        }
        return drugs;
    }

    @NotNull
    private static Drug resolveDrug(@NotNull CkbJsonDatabase ckbJsonDatabase, @NotNull DrugInfo drugInfo) {
        for (JsonDrug drug : ckbJsonDatabase.drugs()) {
            if (drug.id() == drugInfo.id()) {
                return ImmutableDrug.builder()
                        .id(drug.id())
                        .createDate(drug.createDate())
                        .drugName(drug.drugName())
                        .drugClasses(extractDrugClasses(ckbJsonDatabase, drug.drugClasses()))
                        .terms(drug.terms())
                        .synonyms(drug.synonyms())
                        .tradeName(drug.tradeName())
                        .descriptions(extractDrugDescriptions(ckbJsonDatabase, drug.descriptions()))
                        .casRegistryNum(drug.casRegistryNum())
                        .ncitId(drug.ncitId())
                        .build();
            }
        }

        throw new IllegalStateException("Could not resolve CKB drug with id '" + drugInfo.id() + "'");
    }

    @NotNull
    private static List<DrugClass> extractDrugClasses(@NotNull CkbJsonDatabase ckbJsonDatabase,
            @NotNull List<DrugClassInfo> drugClassInfos) {
        List<DrugClass> drugClasses = Lists.newArrayList();
        for (DrugClassInfo drugClassInfo : drugClassInfos) {
            drugClasses.add(resolveDrugClass(ckbJsonDatabase, drugClassInfo));
        }
        return drugClasses;
    }

    @NotNull
    private static DrugClass resolveDrugClass(@NotNull CkbJsonDatabase ckbJsonDatabase, @NotNull DrugClassInfo drugClassInfo) {
        for (JsonDrugClass drugClass : ckbJsonDatabase.drugClasses()) {
            if (drugClassInfo.id() == drugClass.id()) {
                return ImmutableDrugClass.builder()
                        .id(drugClass.id())
                        .createDate(drugClass.createDate())
                        .drugClass(drugClass.drugClass())
                        .build();
            }
        }

        throw new IllegalStateException("Could not resolve CKB drug class with id '" + drugClassInfo.id() + "'");
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