package com.hartwig.hmftools.ckb.interpretation.common.druginterpretation;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodelinterpretation.drug.DrugDescription;
import com.hartwig.hmftools.ckb.datamodelinterpretation.drug.ImmutableDrug;
import com.hartwig.hmftools.ckb.datamodelinterpretation.drug.ImmutableDrugDescription;
import com.hartwig.hmftools.ckb.datamodelinterpretation.drugclass.ImmutableDrugClass;
import com.hartwig.hmftools.ckb.interpretation.common.CommonInterpretationFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.DrugClassInfo;
import com.hartwig.hmftools.ckb.json.common.DrugInfo;
import com.hartwig.hmftools.ckb.json.drug.Drug;
import com.hartwig.hmftools.ckb.json.drugclass.DrugClass;

import org.jetbrains.annotations.NotNull;

public class DrugInterpretationFactory {

    private DrugInterpretationFactory() {

    }

    @NotNull
    public static DrugInterpretation extractDrugInterpretations(@NotNull List<DrugInfo> drugs, @NotNull CkbJsonDatabase ckbEntry) {
        ImmutableDrugInterpretation.Builder outputBuilderDrugInterpretation = ImmutableDrugInterpretation.builder();
        for (DrugInfo drugInfo : drugs) {
            for (Drug drug : ckbEntry.drugs()) {
                if (drugInfo.id() == drug.id()) {
                    outputBuilderDrugInterpretation.addDrugs(ImmutableDrug.builder()
                            .id(drug.id())
                            .drugName(drug.drugName())
                            .terms(drug.terms())
                            .synonyms(drug.synonyms())
                            .tradeName(drug.tradeName())
                            .drugDescriptions(createDrugDescriptions(drug.descriptions(), ckbEntry))
                            .casRegistryNum(drug.casRegistryNum())
                            .ncitId(drug.ncitId())
                            .createDate(drug.createDate())
                            .drugClasses(extractDrugClasses(drug, ckbEntry))
                            .build());
                }
            }
        }
        return outputBuilderDrugInterpretation.build();
    }

    @NotNull
    private static List<com.hartwig.hmftools.ckb.datamodelinterpretation.drugclass.DrugClass> extractDrugClasses(@NotNull Drug drug,
            @NotNull CkbJsonDatabase ckbEntry) {
        List<com.hartwig.hmftools.ckb.datamodelinterpretation.drugclass.DrugClass> drugClasses = Lists.newArrayList();
        for (DrugClassInfo drugClassInfo : drug.drugClasses()) {
            for (DrugClass drugClass : ckbEntry.drugClasses()) {
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
    private static List<DrugDescription> createDrugDescriptions(@NotNull List<DescriptionInfo> descriptionInfos,
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