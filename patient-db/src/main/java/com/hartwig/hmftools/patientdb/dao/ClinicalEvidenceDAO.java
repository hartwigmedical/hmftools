package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CLINICALEVIDENCE;

import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.actionability.EvidenceItem;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep10;

class ClinicalEvidenceDAO {

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
                    CLINICALEVIDENCE.EVENTTYPE,
                    CLINICALEVIDENCE.EVENTMATCH,
                    CLINICALEVIDENCE.NAMEEVIDENCE,
                    CLINICALEVIDENCE.TYPEEVIDENCE,
                    CLINICALEVIDENCE.RESPONSE,
                    CLINICALEVIDENCE.LEVELEVIDENCE,
                    CLINICALEVIDENCE.SOURCEEVIDENCE,
                    CLINICALEVIDENCE.CANCERTYPE,
                    CLINICALEVIDENCE.ISONLABEL);
            items.forEach(trial -> addValues(sample, trial, inserter));
            inserter.execute();
        }
    }

    private static void addValues(@NotNull String sample, @NotNull EvidenceItem evidenceItem, @NotNull InsertValuesStep10 inserter) {
        //noinspection unchecked
        inserter.values(sample,
                evidenceItem.event(),
                evidenceItem.scope().readableString(),
                evidenceItem.drug(),
                evidenceItem.drugsType(),
                evidenceItem.response(),
                evidenceItem.level().readableString(),
                evidenceItem.source().sourceName(),
                evidenceItem.cancerType(),
                evidenceItem.isOnLabel());
    }

    void deleteClinicalEvidenceForSample(@NotNull String sample) {
        context.delete(CLINICALEVIDENCE).where(CLINICALEVIDENCE.SAMPLEID.eq(sample)).execute();
    }
}
