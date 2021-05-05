package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Flagstat.FLAGSTAT;

import com.hartwig.hmftools.common.flagstat.Flagstat;
import com.hartwig.hmftools.common.flagstat.FlagstatQC;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

public class FlagstatDAO {

    @NotNull
    private final DSLContext context;

    FlagstatDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeFlagstats(@NotNull String sample, @NotNull Flagstat refFlagstat, @NotNull Flagstat tumorFlagstat) {
        deleteFlagstatsForSample(sample);

        boolean passQC = FlagstatQC.pass(refFlagstat) && FlagstatQC.pass(tumorFlagstat);
        context.insertInto(FLAGSTAT,
                FLAGSTAT.SAMPLEID,
                FLAGSTAT.REFMAPPEDPROPORTION,
                FLAGSTAT.REFDUPLICATEPROPORTION,
                FLAGSTAT.TUMORMAPPEDPROPORTION,
                FLAGSTAT.TUMORDUPLICATEPROPORTION,
                FLAGSTAT.PASSQC)
                .values(sample,
                        DatabaseUtil.decimal(refFlagstat.mappedProportion()),
                        DatabaseUtil.decimal(refFlagstat.duplicateProportion()),
                        DatabaseUtil.decimal(tumorFlagstat.mappedProportion()),
                        DatabaseUtil.decimal(tumorFlagstat.duplicateProportion()),
                        passQC ? (byte) 1 : (byte) 0)
                .execute();
    }

    void deleteFlagstatsForSample(@NotNull String sample) {
        context.delete(FLAGSTAT).where(FLAGSTAT.SAMPLEID.eq(sample)).execute();
    }
}
