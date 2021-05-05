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
                FLAGSTAT.REFUNIQUEREADCOUNT,
                FLAGSTAT.REFSECONDARYCOUNT,
                FLAGSTAT.REFSUPPLEMENTARYCOUNT,
                FLAGSTAT.REFDUPLICATEPROPORTION,
                FLAGSTAT.REFMAPPEDPROPORTION,
                FLAGSTAT.REFPAIREDINSEQUENCINGPROPORTION,
                FLAGSTAT.REFPROPERLYPAIREDPROPORTION,
                FLAGSTAT.REFWITHITSELFANDMATEMAPPEDPROPORTION,
                FLAGSTAT.REFSINGLETONPROPORTION,
                FLAGSTAT.TUMORUNIQUEREADCOUNT,
                FLAGSTAT.TUMORSECONDARYCOUNT,
                FLAGSTAT.TUMORSUPPLEMENTARYCOUNT,
                FLAGSTAT.TUMORDUPLICATEPROPORTION,
                FLAGSTAT.TUMORMAPPEDPROPORTION,
                FLAGSTAT.TUMORPAIREDINSEQUENCINGPROPORTION,
                FLAGSTAT.TUMORPROPERLYPAIREDPROPORTION,
                FLAGSTAT.TUMORWITHITSELFANDMATEMAPPEDPROPORTION,
                FLAGSTAT.TUMORSINGLETONPROPORTION,
                FLAGSTAT.PASSQC)
                .values(sample,
                        refFlagstat.uniqueReadCount(),
                        refFlagstat.secondaryCount(),
                        refFlagstat.supplementaryCount(),
                        DatabaseUtil.decimal(refFlagstat.duplicateProportion()),
                        DatabaseUtil.decimal(refFlagstat.mappedProportion()),
                        DatabaseUtil.decimal(refFlagstat.pairedInSequencingProportion()),
                        DatabaseUtil.decimal(refFlagstat.properlyPairedProportion()),
                        DatabaseUtil.decimal(refFlagstat.withItselfAndMateMappedProportion()),
                        DatabaseUtil.decimal(refFlagstat.singletonProportion()),
                        tumorFlagstat.uniqueReadCount(),
                        tumorFlagstat.secondaryCount(),
                        tumorFlagstat.supplementaryCount(),
                        DatabaseUtil.decimal(tumorFlagstat.duplicateProportion()),
                        DatabaseUtil.decimal(tumorFlagstat.mappedProportion()),
                        DatabaseUtil.decimal(tumorFlagstat.pairedInSequencingProportion()),
                        DatabaseUtil.decimal(tumorFlagstat.properlyPairedProportion()),
                        DatabaseUtil.decimal(tumorFlagstat.withItselfAndMateMappedProportion()),
                        DatabaseUtil.decimal(tumorFlagstat.singletonProportion()),
                        passQC ? (byte) 1 : (byte) 0)
                .execute();
    }

    void deleteFlagstatsForSample(@NotNull String sample) {
        context.delete(FLAGSTAT).where(FLAGSTAT.SAMPLEID.eq(sample)).execute();
    }
}
