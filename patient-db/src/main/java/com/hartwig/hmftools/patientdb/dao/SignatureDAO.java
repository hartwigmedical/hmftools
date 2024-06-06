package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SIGNATURE;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sigs.ImmutableSignatureAllocation;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep5;
import org.jooq.Record;
import org.jooq.Result;

class SignatureDAO
{
    private final DSLContext context;

    SignatureDAO(final DSLContext context)
    {
        this.context = context;
    }

    public void write(final String sample, final List<SignatureAllocation> sigAllocations)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(SIGNATURE).where(SIGNATURE.SAMPLEID.eq(sample)).execute();

        for (List<SignatureAllocation> batch : Iterables.partition(sigAllocations, DB_BATCH_INSERT_SIZE))
        {
            InsertValuesStep5 inserter = context.insertInto(SIGNATURE,
                    SIGNATURE.SAMPLEID,
                    SIGNATURE.MODIFIED,
                    SIGNATURE.SIGNATURE_,
                    SIGNATURE.ALLOCATION,
                    SIGNATURE.PERCENT);

            batch.forEach(entry -> addRecord(timestamp, inserter, sample, entry));
            inserter.execute();
        }
    }

    private static void addRecord(final Timestamp timestamp, final InsertValuesStep5 inserter, final String sample,
            final SignatureAllocation sigAllocation)
    {
        inserter.values(sample,
                timestamp,
                sigAllocation.signature(),
                DatabaseUtil.decimal(sigAllocation.allocation()),
                DatabaseUtil.decimal(sigAllocation.percent()));
    }

    public List<SignatureAllocation> readAllocations(final String sample)
    {
        List<SignatureAllocation> sigAllocationList = Lists.newArrayList();

        Result<Record> result = context.select().from(SIGNATURE).where(SIGNATURE.SAMPLEID.eq(sample)).fetch();

        for (Record record : result)
        {
            SignatureAllocation sigAllocation = ImmutableSignatureAllocation.builder()
                    .signature(record.getValue(SIGNATURE.SIGNATURE_))
                    .allocation(record.getValue(SIGNATURE.ALLOCATION))
                    .percent(record.getValue(SIGNATURE.PERCENT))
                    .build();

            sigAllocationList.add(sigAllocation);
        }

        return sigAllocationList;
    }

    void deleteSignatureDataForSample(final String sample)
    {
        context.delete(SIGNATURE).where(SIGNATURE.SAMPLEID.eq(sample)).execute();
    }
}
