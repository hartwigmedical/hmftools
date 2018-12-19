package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CLINICALEVIDENCE;

import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.actionability.EvidenceItem;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep10;
import org.jooq.InsertValuesStep11;
import org.jooq.InsertValuesStep19;

public class ClinicalEvidenceDAO {

    @NotNull
    private final DSLContext context;

    ClinicalEvidenceDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeClinicalEvidence(@NotNull String sample, @NotNull List<EvidenceItem> evidenceItem) {
        deleteClinicalEvidenceForSample(sample);

        for (List<EvidenceItem> items : Iterables.partition(evidenceItem, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep10 inserter = context.insertInto(CLINICALEVIDENCE,
                    CLINICALEVIDENCE.SAMPLEID,
//                    CLINICALEVIDENCE.TYPEVARIANT,
//                    CLINICALEVIDENCE.GENE,
//                    CLINICALEVIDENCE.CHOMOSOME,
//                    CLINICALEVIDENCE.POSITION,
//                    CLINICALEVIDENCE.REF,
//                    CLINICALEVIDENCE.ALT,
//                    CLINICALEVIDENCE.CNVTYPE,
//                    CLINICALEVIDENCE.FUSIONFIVEGENE,
//                    CLINICALEVIDENCE.FUSIONTHREEGENE,
                    CLINICALEVIDENCE.EVENTTYPE,
                    CLINICALEVIDENCE.EVENTMATCH,
                    CLINICALEVIDENCE.DRUG,
                    CLINICALEVIDENCE.DRUGSTYPE,
                    CLINICALEVIDENCE.RESPONSE,
                    CLINICALEVIDENCE.CANCERTYPE,
                    CLINICALEVIDENCE.LABEL,
                    CLINICALEVIDENCE.EVIDENCELEVEL,
                    CLINICALEVIDENCE.EVIDENCESOURCE);
            items.forEach(trial -> addValues(sample, trial, inserter));
            inserter.execute();
        }

    }

    private static void addValues(@NotNull String sample, @NotNull EvidenceItem evidenceItem, @NotNull InsertValuesStep10 inserter) {
        //noinspection unchecked
        inserter.values(sample,
//                evidenceItem.type(),
//                evidenceItem.gene(),
//                evidenceItem.chromosome(),
//                evidenceItem.position(),
//                evidenceItem.ref(),
//                evidenceItem.alt(),
//                evidenceItem.cnvType(),
//                evidenceItem.fusionFiveGene(),
//                evidenceItem.fusionThreeGene(),
                evidenceItem.event(),
                evidenceItem.scope().readableString(),
                evidenceItem.drug(),
                evidenceItem.drugsType(),
                evidenceItem.response(),
                evidenceItem.cancerType(),
                evidenceItem.isOnLabel() ? "Tumor Type specific" : "Other tumor types specific",
                evidenceItem.level().readableString(),
                evidenceItem.source().sourceName());

    }

    void deleteClinicalEvidenceForSample(@NotNull String sample) {
        context.delete(CLINICALEVIDENCE).where(CLINICALEVIDENCE.SAMPLEID.eq(sample)).execute();
    }
}
