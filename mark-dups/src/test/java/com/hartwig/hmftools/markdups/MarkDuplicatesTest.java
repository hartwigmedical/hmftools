package com.hartwig.hmftools.markdups;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.junit.Test;

public class MarkDuplicatesTest
{
    private final ReadIdGenerator mReadIdGen;
    private final MarkDupsConfig mConfig;
    private final RecordWriter mWriter;
    private final PartitionDataStore mPartitionDataStore;

    public MarkDuplicatesTest()
    {
        mReadIdGen = new ReadIdGenerator();
        mConfig = new MarkDupsConfig();
        mWriter = new RecordWriter(mConfig);
        mPartitionDataStore = new PartitionDataStore(mConfig);

    }

    @Test
    public void testSinglePartition()
    {
        ChrBaseRegion chrBaseRegion = new ChrBaseRegion(CHR_1, 1, 100000);
        ChromosomeReader chromosomeReader = new ChromosomeReader(chrBaseRegion, mConfig, mWriter, mPartitionDataStore);

        // things to test:
        // - paired and unpaired reads
        // - supplementaries and mates coming first
        // -




    }


}
