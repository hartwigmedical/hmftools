package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PROTECT;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.protect.ProtectEvidenceItem;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep10;

class ProtectDAO {

    @NotNull
    private final DSLContext context;

    ProtectDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void write(@NotNull String sample, @NotNull List<ProtectEvidenceItem> evidence) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        deleteEvidenceForSample(sample);
        for (List<ProtectEvidenceItem> batch : Iterables.partition(evidence, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep10 inserter = context.insertInto(PROTECT,
                    PROTECT.SAMPLEID,
                    PROTECT.EVENT,
                    PROTECT.SOURCE,
                    PROTECT.TREATMENT,
                    PROTECT.LEVEL,
                    PROTECT.DIRECTION,
                    PROTECT.URLS,
                    PROTECT.ONLABEL,
                    PROTECT.REPORTED,
                    PROTECT.MODIFIED);
            batch.forEach(entry -> addRecord(timestamp, inserter, sample, entry));
            inserter.execute();
        }
    }

    private static void addRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStep10 inserter, @NotNull String sample,
            @NotNull ProtectEvidenceItem evidence) {
        StringJoiner urlJoiner = new StringJoiner(",");
        for (String url : evidence.urls()) {
            urlJoiner.add(url);
        }

        inserter.values(sample,
                evidence.genomicEvent(),
                evidence.source().toString(),
                evidence.treatment(),
                evidence.level().toString(),
                evidence.direction().toString(),
                urlJoiner.toString(),
                evidence.onLabel(),
                evidence.reported(),
                timestamp);
    }

    void deleteEvidenceForSample(@NotNull String sample) {
        context.delete(PROTECT).where(PROTECT.SAMPLEID.eq(sample)).execute();
    }
}
