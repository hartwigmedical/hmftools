package com.hartwig.hmftools.patientdb;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.common.metrics.ImmutableBamFlagStats;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.FlagstatRecord;

import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class LoadFlagstatDataTest extends DatabaseTestBase
{
    @Test
    public void canWriteBamFlagStats()
    {
        BamFlagStats flagStatsTemplate = ImmutableBamFlagStats.builder()
                .uniqueReadCount(0)
                .secondaryCount(0)
                .supplementaryCount(0)
                .duplicateProportion(0)
                .mappedProportion(0)
                .pairedInSequencingProportion(0)
                .properlyPairedProportion(0)
                .withItselfAndMateMappedProportion(0)
                .singletonProportion(0)
                .build();

        BamFlagStats refFlagStats = flagStatsTemplate;
        BamFlagStats tumorFlagStats = flagStatsTemplate;

        databaseAccess.writeFlagstats(TEST_SAMPLE_ID, refFlagStats, tumorFlagStats);

        List<FlagstatRecord> flagStatRecords = fetchTable(Tables.FLAGSTAT, FlagstatRecord.class);
        assertEquals(1, flagStatRecords.size());
    }
}
