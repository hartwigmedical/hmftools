package com.hartwig.hmftools.bamtools.btofmc.readcache;

import static com.hartwig.hmftools.bamtools.btofmc.TestUtil.mappedPrimarySamRecord;
import static com.hartwig.hmftools.bamtools.btofmc.TestUtil.pairedMappedPrimarySamRecords;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.function.BiConsumer;

import com.google.inject.AbstractModule;
import com.google.inject.Guice;
import com.google.inject.Injector;
import com.google.inject.Singleton;
import com.google.inject.Stage;
import com.google.inject.TypeLiteral;
import com.hartwig.hmftools.bamtools.btofmc.BamToFastqConfig;
import com.hartwig.hmftools.bamtools.btofmc.TestUtil.ReadPairCollector;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ReadCacheTest
{
    private static class TestModule extends AbstractModule
    {
        private final BamToFastqConfig mConfig;

        public TestModule()
        {
            mConfig = new BamToFastqConfig(false);
        }

        @Override
        public void configure()
        {
            bind(BamToFastqConfig.class).toInstance(mConfig);

            bind(ReadCache.class).in(Singleton.class);
            bind(ReadCacheInterface.class).to(ReadCache.class).in(Singleton.class);

            bind(ReadPairCollector.class).in(Singleton.class);
            bind(new TypeLiteral<BiConsumer<SAMRecord, SAMRecord>>()
            {
            }).to(ReadPairCollector.class).in(Singleton.class);
        }
    }

    @Test
    public void testSizeOnEmpty()
    {
        Injector injector = Guice.createInjector(Stage.PRODUCTION, new TestModule());
        ReadCacheInterface readCache = injector.getInstance(ReadCacheInterface.class);
        ReadPairCollector readPairCollector = injector.getInstance(ReadPairCollector.class);

        assertEquals(0, readCache.size());
        assertTrue(readCache.isEmpty());
        assertTrue(readPairCollector.ReadPairs.isEmpty());
    }

    @Test
    public void testSizeOnAfterAddingSingleRead()
    {
        Injector injector = Guice.createInjector(Stage.PRODUCTION, new TestModule());
        ReadCacheInterface readCache = injector.getInstance(ReadCacheInterface.class);
        ReadPairCollector readPairCollector = injector.getInstance(ReadPairCollector.class);

        SAMRecord read = mappedPrimarySamRecord("READ_01", CHR_1, 1);
        readCache.add(read);

        assertEquals(1, readCache.size());
        assertFalse(readCache.isEmpty());
        assertTrue(readPairCollector.ReadPairs.isEmpty());
    }

    @Test
    public void testSizeOnAfterAddingPairedReads()
    {
        Injector injector = Guice.createInjector(Stage.PRODUCTION, new TestModule());
        ReadCacheInterface readCache = injector.getInstance(ReadCacheInterface.class);
        ReadPairCollector readPairCollector = injector.getInstance(ReadPairCollector.class);

        Pair<SAMRecord, SAMRecord> readPair = pairedMappedPrimarySamRecords("READ_01", CHR_1, 1, CHR_2, 1);
        readCache.add(readPair.getLeft());
        readCache.add(readPair.getRight());

        assertEquals(0, readCache.size());
        assertTrue(readCache.isEmpty());
        assertEquals(1, readPairCollector.ReadPairs.size());
    }
}
