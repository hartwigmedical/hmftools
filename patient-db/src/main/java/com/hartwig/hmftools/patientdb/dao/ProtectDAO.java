package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PROTECT;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectSource;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep19;

class ProtectDAO {

    @NotNull
    private final DSLContext context;

    ProtectDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void write(@NotNull String sample, @NotNull List<ProtectEvidence> evidence) {
        deleteEvidenceForSample(sample);

        Timestamp timestamp = new Timestamp(new Date().getTime());
        for (List<ProtectEvidence> batch : Iterables.partition(evidence, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep19 inserter = context.insertInto(PROTECT,
                    PROTECT.SAMPLEID,
                    PROTECT.GENE,
                    PROTECT.TRANSCRIPT,
                    PROTECT.ISCANONICAL,
                    PROTECT.EVENT,
                    PROTECT.EVENTISHIGHDRIVER,
                    PROTECT.GERMLINE,
                    PROTECT.REPORTED,
                    PROTECT.TREATMENT,
                    PROTECT.ONLABEL,
                    PROTECT.LEVEL,
                    PROTECT.DIRECTION,
                    PROTECT.SOURCE,
                    PROTECT.SOURCEEVENT,
                    PROTECT.SOURCEURLS,
                    PROTECT.EVIDENCETYPE,
                    PROTECT.RANGERANK,
                    PROTECT.EVIDENCEURLS,
                    PROTECT.MODIFIED);
            batch.forEach(entry -> addRecord(timestamp, inserter, sample, entry));
            inserter.execute();
        }
    }

    private static void addRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStep19 inserter, @NotNull String sample,
            @NotNull ProtectEvidence evidence) {
        StringJoiner sourceUrlJoiner = new StringJoiner(",");
        StringJoiner evidenceUrlJoiner = new StringJoiner(",");

        for (ProtectSource protectSource : evidence.protectSources()) {
            for (String sourceUrl : protectSource.sourceUrls()) {
                sourceUrlJoiner.add(sourceUrl);
            }

            for (String evidenceUrl : protectSource.evidenceUrls()) {
                evidenceUrlJoiner.add(evidenceUrl);
            }
            String knowledgebase = protectSource.source().technicalDisplay();
            String evidenceType = protectSource.evidenceType().display();
            Integer rank = protectSource.rangeRank();

            inserter.values(sample,
                    evidence.gene(),
                    evidence.transcript(),
                    evidence.isCanonical(),
                    evidence.event(),
                    evidence.eventIsHighDriver(),
                    evidence.germline(),
                    evidence.reported(),
                    evidence.treatment(),
                    evidence.onLabel(),
                    evidence.level().toString(),
                    evidence.direction().toString(),
                    knowledgebase,
                    protectSource.sourceEvent(),
                    sourceUrlJoiner.toString(),
                    evidenceType,
                    rank,
                    evidenceUrlJoiner.toString(),
                    timestamp);
        }
    }

    void deleteEvidenceForSample(@NotNull String sample) {
        context.delete(PROTECT).where(PROTECT.SAMPLEID.eq(sample)).execute();
    }
}