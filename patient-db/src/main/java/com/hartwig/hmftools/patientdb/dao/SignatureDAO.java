package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SIGNATURE;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.Field;
import org.jooq.InsertValuesStep5;

class SignatureDAO
{
    @NotNull
    private final DSLContext context;

    SignatureDAO(@NotNull final DSLContext context)
    {
        this.context = context;
    }

    public void write(@NotNull String sample, @NotNull List<SignatureAllocation> sigAllocations)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(SIGNATURE).where(SIGNATURE.SAMPLEID.eq(sample)).execute();

        for (List<SignatureAllocation> batch : Iterables.partition(sigAllocations, DB_BATCH_INSERT_SIZE))
        {
            InsertValuesStep5 inserter = context.insertInto(SIGNATURE,
                    SIGNATURE.SAMPLEID,
                    SIGNATURE.MODIFIED,
                    (Field<?>) SIGNATURE.SIGNATURE,
                    SIGNATURE.ALLOCATION,
                    SIGNATURE.PERCENT);

            batch.forEach(entry -> addRecord(timestamp, inserter, sample, entry));
            inserter.execute();
        }
    }

    private static void addRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStep5 inserter, @NotNull String sample,
            @NotNull SignatureAllocation sigAllocation)
    {
        inserter.values(sample,
                timestamp,
                sigAllocation.signature(),
                sigAllocation.allocation(),
                sigAllocation.percent());
    }

    void deleteSignatureDataForSample(@NotNull String sample)
    {
        context.delete(SIGNATURE).where(SIGNATURE.SAMPLEID.eq(sample)).execute();
    }
}
