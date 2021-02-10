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

public class DrugsInterpretationFactory {

    private DrugsInterpretationFactory() {

    }

    @NotNull
    public static DrugsInterpretation extractDrugsInterpretation(@NotNull List<DrugInfo> drugs, @NotNull CkbJsonDatabase ckbEntry) {
        ImmutableDrugsInterpretation.Builder outputBuilderDrugInterpretation = ImmutableDrugsInterpretation.builder();
        for (DrugInfo drugInfo : drugs) {
            for (Drug drug : ckbEntry.drugs()) {
                if (drugInfo.id() == drug.id()) {
                    outputBuilderDrugInterpretation.drug(ImmutableDrug.builder()
                            .id(drug.id())
                            .drugName(drug.drugName())
                            .terms(drug.terms())
                            .synonyms(drug.synonyms())
                            .tradeName(drug.tradeName())
                            .drugDescriptions(createDrugDescriptions(drug.descriptions(), ckbEntry))
                            .casRegistryNum(drug.casRegistryNum())
                            .ncitId(drug.ncitId())
                            .createDate(drug.createDate())
                            .build());

                    for (DrugClassInfo drugClassInfo : drug.drugClasses()) {
                        for (DrugClass drugClass : ckbEntry.drugClasses()) {
                            if (drugClassInfo.id() == drugClass.id()) {
                                outputBuilderDrugInterpretation.drugClass(ImmutableDrugClass.builder()
                                        .id(drugClass.id())
                                        .drugClass(drugClass.drugClass())
                                        .createDate(drugClass.createDate())
                                        .build());
                            }
                        }
                    }
                }
            }
        }
        return outputBuilderDrugInterpretation.build();
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