package com.hartwig.hmftools.serve.dao;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep7;

import static com.hartwig.hmftools.serve.database.tables.Actin.ACTIN;

public class ActinDAO {

    @NotNull
    private final DSLContext context;

    ActinDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void write(@NotNull List<ActinEntry> trials) {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        for (List<ActinEntry> batch : Iterables.partition(trials, Utils.DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep7 inserter = context.insertInto(ACTIN,
                    ACTIN.MODIFIED,
                    ACTIN.TRIAL,
                    ACTIN.COHORT,
                    ACTIN.RULE,
                    ACTIN.GENE,
                    ACTIN.MUTATION,
                    ACTIN.ISUSEDASINCLUSION);
            batch.forEach(entry -> addRecord(timestamp, inserter, entry));
            inserter.execute();
        }

    }

    private static void addRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStep7 inserter,
            @NotNull ActinEntry trial) {
        inserter.values(timestamp,
                trial.trial(),
                trial.cohort(),
                trial.rule(),
                trial.gene(),
                trial.mutation(),
                trial.isUsedAsInclusion());
    }
}
