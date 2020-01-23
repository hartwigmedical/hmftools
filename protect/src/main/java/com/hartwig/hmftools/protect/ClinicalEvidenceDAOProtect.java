package com.hartwig.hmftools.protect;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CLINICALEVIDENCEPROTECT;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.protect.actionability.EvidenceItem;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep10;

public class ClinicalEvidenceDAOProtect {

    @NotNull
    private final DSLContext context;

    ClinicalEvidenceDAOProtect(@NotNull final DSLContext context) {
        this.context = context;
    }

    public void writeClinicalEvidence(@NotNull String sample, @NotNull List<EvidenceItem> evidenceItem) {
        deleteClinicalEvidenceForSample(sample);

        Timestamp timestamp = new Timestamp(new Date().getTime());

        for (List<EvidenceItem> items : Iterables.partition(evidenceItem, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep10 inserter = context.insertInto(CLINICALEVIDENCEPROTECT,
                    CLINICALEVIDENCEPROTECT.SAMPLEID,
                    CLINICALEVIDENCEPROTECT.MODIFIED,
                    CLINICALEVIDENCEPROTECT.EVENT,
                    CLINICALEVIDENCEPROTECT.NAME,
                    CLINICALEVIDENCEPROTECT.TYPE,
                    CLINICALEVIDENCEPROTECT.RESPONSE,
                    CLINICALEVIDENCEPROTECT.LEVEL,
                    CLINICALEVIDENCEPROTECT.SOURCE,
                    CLINICALEVIDENCEPROTECT.CANCERTYPE,
                    CLINICALEVIDENCEPROTECT.ISONLABEL);
            items.forEach(trial -> addValues(sample, trial, inserter, timestamp));
            inserter.execute();
        }
    }

    private static void addValues(@NotNull String sample, @NotNull EvidenceItem evidenceItem, @NotNull InsertValuesStep10 inserter,
            @NotNull Timestamp timestamp) {
        inserter.values(sample,
                timestamp,
                evidenceItem.event(),
                evidenceItem.drug(),
                evidenceItem.drugsType(),
                evidenceItem.response(),
                evidenceItem.level().readableString(),
                evidenceItem.source().sourceName(),
                evidenceItem.cancerType(),
                evidenceItem.isOnLabel());
    }

    public void deleteClinicalEvidenceForSample(@NotNull String sample) {
        context.delete(CLINICALEVIDENCEPROTECT).where(CLINICALEVIDENCEPROTECT.SAMPLEID.eq(sample)).execute();
    }
}
