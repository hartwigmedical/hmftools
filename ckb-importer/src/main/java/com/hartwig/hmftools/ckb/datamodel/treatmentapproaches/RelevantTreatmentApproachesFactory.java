package com.hartwig.hmftools.ckb.datamodel.treatmentapproaches;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.drug.DrugFactory;
import com.hartwig.hmftools.ckb.datamodel.reference.ReferenceFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.TreatmentApproachInfo;
import com.hartwig.hmftools.ckb.json.treatmentapproach.JsonTreatmentApproach;

import org.jetbrains.annotations.NotNull;

public final class RelevantTreatmentApproachesFactory {

    private RelevantTreatmentApproachesFactory() {
    }

    public static List<RelevantTreatmentApproaches> extractRelevantTreatmentApproaches(@NotNull CkbJsonDatabase ckbJsonDatabase,
            @NotNull List<TreatmentApproachInfo> treatmentApproachInfos) {
        List<RelevantTreatmentApproaches> relevantTreatmentApproach = Lists.newArrayList();

        for (TreatmentApproachInfo treatmentApproachInfo : treatmentApproachInfos) {
            RelevantTreatmentApproaches resolvedRelevantTreatmentApproaches =
                    resolveRelevantTreatmentApproaches(ckbJsonDatabase, treatmentApproachInfo);
            relevantTreatmentApproach.add(resolvedRelevantTreatmentApproaches);
        }
        return relevantTreatmentApproach;
    }

    @NotNull
    private static RelevantTreatmentApproaches resolveRelevantTreatmentApproaches(@NotNull CkbJsonDatabase ckbJsonDatabase,
            @NotNull TreatmentApproachInfo treatmentApproachInfo) {
        for (JsonTreatmentApproach treatmentApproah : ckbJsonDatabase.treatmentApproaches()) {
            if (treatmentApproah.id() == treatmentApproachInfo.id()) {
                return ImmutableRelevantTreatmentApproaches.builder()
                        .id(treatmentApproah.id())
                        .drugClass(treatmentApproah.drugClass() != null ? DrugFactory.resolveDrugClass(ckbJsonDatabase,
                                treatmentApproah.drugClass()) : null)
                        .references(ReferenceFactory.extractReferences(ckbJsonDatabase, treatmentApproah.references()))
                        .createDate(treatmentApproah.createDate())
                        .updateDate(treatmentApproah.updateDate())
                        .build();
            }
        }
        throw new IllegalStateException("Could not resolve CKB treatment approach with id '" + treatmentApproachInfo.id() + "'");
    }
}