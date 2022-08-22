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
import org.jooq.InsertValuesStep6;
import org.jooq.Record;
import org.jooq.Result;

class SignatureDAO {
    @NotNull
    private final DSLContext context;

    SignatureDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    public void write(@NotNull String sample, @NotNull String isolationBarcode, @NotNull List<SignatureAllocation> sigAllocations) {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(SIGNATURE).where(SIGNATURE.SAMPLEID.eq(sample)).execute();

        for (List<SignatureAllocation> batch : Iterables.partition(sigAllocations, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep6 inserter = context.insertInto(SIGNATURE,
                    SIGNATURE.SAMPLEID,
                    SIGNATURE.ISOLATIONBARCODE,
                    SIGNATURE.MODIFIED,
                    SIGNATURE.SIGNATURE_,
                    SIGNATURE.ALLOCATION,
                    SIGNATURE.PERCENT);

            batch.forEach(entry -> addRecord(timestamp, inserter, sample, isolationBarcode, entry));
            inserter.execute();
        }
    }

    private static void addRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStep6 inserter, @NotNull String sample,
            @NotNull String isolationBarcode, @NotNull SignatureAllocation sigAllocation) {
        inserter.values(sample,
                isolationBarcode,
                timestamp,
                sigAllocation.signature(),
                DatabaseUtil.decimal(sigAllocation.allocation()),
                DatabaseUtil.decimal(sigAllocation.percent()));
    }

    @NotNull
    public List<SignatureAllocation> readAllocations(@NotNull String sample) {
        List<SignatureAllocation> sigAllocationList = Lists.newArrayList();

        Result<Record> result = context.select().from(SIGNATURE).where(SIGNATURE.SAMPLEID.eq(sample)).fetch();

        for (Record record : result) {
            SignatureAllocation sigAllocation = ImmutableSignatureAllocation.builder()
                    .signature(record.getValue(SIGNATURE.SIGNATURE_))
                    .allocation(record.getValue(SIGNATURE.ALLOCATION))
                    .percent(record.getValue(SIGNATURE.PERCENT))
                    .build();

            sigAllocationList.add(sigAllocation);
        }

        return sigAllocationList;
    }

    void deleteSignatureDataForSample(@NotNull String sample) {
        context.delete(SIGNATURE).where(SIGNATURE.SAMPLEID.eq(sample)).execute();
    }
}
