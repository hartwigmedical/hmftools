package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.BQR;

import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.patientdb.bqr.BQREntry;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep8;

public class BQRDAO {

    @NotNull
    private final DSLContext context;

    BQRDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void write(@NotNull String sample, @NotNull List<BQREntry> refEntries, @NotNull List<BQREntry> tumorEntries) {
        deleteBQRForSample(sample);

        InsertValuesStep8 inserter = context.insertInto(BQR,
                BQR.SAMPLEID,
                BQR.SAMPLETYPE,
                BQR.REF,
                BQR.ALT,
                BQR.TRINUCLEOTIDECONTEXT,
                BQR.BQRCOUNT,
                BQR.ORIGQUALITY,
                BQR.RECALIBRATEDQUALITY);

        for (List<BQREntry> batch : Iterables.partition(refEntries, DB_BATCH_INSERT_SIZE)) {
            batch.forEach(entry -> addRecord(inserter, sample, "REF", entry));
            inserter.execute();
        }

        for (List<BQREntry> batch : Iterables.partition(tumorEntries, DB_BATCH_INSERT_SIZE)) {
            batch.forEach(entry -> addRecord(inserter, sample, "TUMOR", entry));
            inserter.execute();
        }
    }

    private static void addRecord(@NotNull InsertValuesStep8 inserter, @NotNull String sample, @NotNull String sampleType,
            @NotNull BQREntry entry) {
        inserter.values(sample,
                sampleType,
                entry.ref(),
                entry.alt(),
                entry.trinucleotideContext(),
                entry.count(),
                entry.origQuality(),
                entry.recalibratedQuality());
    }

    void deleteBQRForSample(@NotNull String sample) {
        context.delete(BQR).where(BQR.SAMPLEID.eq(sample)).execute();
    }
}
