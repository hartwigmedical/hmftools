package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Copynumber.COPYNUMBER;

import java.sql.Timestamp;
import java.util.Collections;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep9;
import org.jooq.Record;
import org.jooq.Result;

class CopyNumberDAO {

    @NotNull
    private final DSLContext context;

    CopyNumberDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    @NotNull
    public List<PurpleCopyNumber> read(@NotNull final String sample) {
        List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();

        Result<Record> result = context.select().from(COPYNUMBER).where(COPYNUMBER.SAMPLEID.eq(sample)).fetch();

        for (Record record : result) {
            copyNumbers.add(ImmutablePurpleCopyNumber.builder()
                    .chromosome(record.getValue(COPYNUMBER.CHROMOSOME))
                    .start(record.getValue(COPYNUMBER.START))
                    .end(record.getValue(COPYNUMBER.END))
                    .bafCount(record.getValue(COPYNUMBER.BAFCOUNT))
                    .averageActualBAF(record.getValue(COPYNUMBER.ACTUALBAF))
                    .averageObservedBAF(record.getValue(COPYNUMBER.OBSERVEDBAF))
                    .averageTumorCopyNumber(record.getValue(COPYNUMBER.COPYNUMBER_))
                    .build());
        }

        Collections.sort(copyNumbers);
        return copyNumbers;
    }

    public void write(@NotNull final String sample, @NotNull List<PurpleCopyNumber> copyNumbers) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(COPYNUMBER).where(COPYNUMBER.SAMPLEID.eq(sample)).execute();

        InsertValuesStep9 inserter = context.insertInto(COPYNUMBER, COPYNUMBER.SAMPLEID, COPYNUMBER.CHROMOSOME,
                COPYNUMBER.START, COPYNUMBER.END, COPYNUMBER.BAFCOUNT, COPYNUMBER.OBSERVEDBAF, COPYNUMBER.ACTUALBAF,
                COPYNUMBER.COPYNUMBER_, COPYNUMBER.MODIFIED);

        copyNumbers.forEach(x -> addCopynumberRecord(timestamp, inserter, sample, x));
        inserter.execute();
    }

    private void addCopynumberRecord(Timestamp timestamp, InsertValuesStep9 inserter, String sample,
            PurpleCopyNumber region) {
        inserter.values(sample, region.chromosome(), region.start(), region.end(), region.bafCount(),
                region.averageObservedBAF(), region.averageActualBAF(), region.averageTumorCopyNumber(), timestamp);
    }
}
