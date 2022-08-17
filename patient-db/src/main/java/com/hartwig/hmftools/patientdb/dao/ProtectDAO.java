package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PROTECT;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.KnowledgebaseSource;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep19;
import org.jooq.InsertValuesStep21;

class ProtectDAO {

    private static final String TREATMENT_APPROACH_DELIMITER = ",";

    @NotNull
    private final DSLContext context;

    ProtectDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void write(@NotNull String sample, @NotNull List<ProtectEvidence> evidence) {
        deleteEvidenceForSample(sample);

        Timestamp timestamp = new Timestamp(new Date().getTime());
        for (List<ProtectEvidence> batch : Iterables.partition(evidence, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep21 inserter = context.insertInto(PROTECT,
                    PROTECT.SAMPLEID,
                    PROTECT.GENE,
                    PROTECT.TRANSCRIPT,
                    PROTECT.ISCANONICAL,
                    PROTECT.EVENT,
                    PROTECT.EVENTISHIGHDRIVER,
                    PROTECT.GERMLINE,
                    PROTECT.REPORTED,
                    PROTECT.TREATMENT,
                    PROTECT.SOURCETREATMENTAPPROACH,
                    PROTECT.TREATMENTAPPROACH,
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

    private static void addRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStep21 inserter, @NotNull String sample,
            @NotNull ProtectEvidence evidence) {
        for (KnowledgebaseSource source : evidence.sources()) {
            StringJoiner sourceUrlJoiner = new StringJoiner(",");
            for (String sourceUrl : source.sourceUrls()) {
                sourceUrlJoiner.add(sourceUrl);
            }

            StringJoiner evidenceUrlJoiner = new StringJoiner(",");
            for (String evidenceUrl : source.evidenceUrls()) {
                evidenceUrlJoiner.add(evidenceUrl);
            }

            inserter.values(sample,
                    evidence.gene(),
                    evidence.transcript(),
                    evidence.isCanonical(),
                    evidence.event(),
                    evidence.eventIsHighDriver(),
                    evidence.germline(),
                    evidence.reported(),
                    evidence.treatment().treament(),
                    treatmentApproachToString(evidence.treatment().sourceRelevantTreatmentApproaches()),
                    treatmentApproachToString(evidence.treatment().relevantTreatmentApproaches()),
                    evidence.onLabel(),
                    evidence.level().toString(),
                    evidence.direction().toString(),
                    source.name().toString(),
                    source.sourceEvent(),
                    sourceUrlJoiner.toString().isEmpty() ? null : sourceUrlJoiner.toString(),
                    source.evidenceType().toString(),
                    source.rangeRank(),
                    evidenceUrlJoiner.toString().isEmpty() ? null : evidenceUrlJoiner.toString(),
                    timestamp);
        }
    }

    @Nullable
    public static String treatmentApproachToString(@NotNull Set<String> treatmentApproaches) {
        StringJoiner joiner = new StringJoiner(TREATMENT_APPROACH_DELIMITER);
        for (String url : treatmentApproaches) {
            joiner.add(url);
        }
        return joiner.toString().isEmpty() ? null : joiner.toString();
    }

    void deleteEvidenceForSample(@NotNull String sample) {
        context.delete(PROTECT).where(PROTECT.SAMPLEID.eq(sample)).execute();
    }
}