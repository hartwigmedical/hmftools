package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.common.metrics.WGSMetricQC.hasSufficientMappedProportion;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Flagstat.FLAGSTAT;

import com.hartwig.hmftools.common.metrics.BamFlagStats;

import org.jooq.DSLContext;

public class FlagstatDAO
{
    private final DSLContext mContext;

    public FlagstatDAO(final DSLContext context)
    {
        this.mContext = context;
    }

    public void writeFlagstats(final String sample, final BamFlagStats refFlagStats, final BamFlagStats tumorFlagStats)
    {
        deleteFlagstatsForSample(sample);

        boolean passQC = hasSufficientMappedProportion(refFlagStats) && hasSufficientMappedProportion(tumorFlagStats);
        
        mContext.insertInto(FLAGSTAT,
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
                        refFlagStats.uniqueReadCount(),
                        refFlagStats.secondaryCount(),
                        refFlagStats.supplementaryCount(),
                        DatabaseUtil.decimal(refFlagStats.duplicateProportion()),
                        DatabaseUtil.decimal(refFlagStats.mappedProportion()),
                        DatabaseUtil.decimal(refFlagStats.pairedInSequencingProportion()),
                        DatabaseUtil.decimal(refFlagStats.properlyPairedProportion()),
                        DatabaseUtil.decimal(refFlagStats.withItselfAndMateMappedProportion()),
                        DatabaseUtil.decimal(refFlagStats.singletonProportion()),
                        tumorFlagStats.uniqueReadCount(),
                        tumorFlagStats.secondaryCount(),
                        tumorFlagStats.supplementaryCount(),
                        DatabaseUtil.decimal(tumorFlagStats.duplicateProportion()),
                        DatabaseUtil.decimal(tumorFlagStats.mappedProportion()),
                        DatabaseUtil.decimal(tumorFlagStats.pairedInSequencingProportion()),
                        DatabaseUtil.decimal(tumorFlagStats.properlyPairedProportion()),
                        DatabaseUtil.decimal(tumorFlagStats.withItselfAndMateMappedProportion()),
                        DatabaseUtil.decimal(tumorFlagStats.singletonProportion()),
                        passQC ? (byte) 1 : (byte) 0)
                .execute();
    }

    public void deleteFlagstatsForSample(final String sample)
    {
        mContext.delete(FLAGSTAT).where(FLAGSTAT.SAMPLEID.eq(sample)).execute();
    }
}
