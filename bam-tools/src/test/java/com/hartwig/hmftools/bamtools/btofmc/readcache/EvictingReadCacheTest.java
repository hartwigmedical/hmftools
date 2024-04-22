package com.hartwig.hmftools.bamtools.btofmc.readcache;

import static com.hartwig.hmftools.bamtools.btofmc.TestUtil.mappedPrimarySamRecord;
import static com.hartwig.hmftools.bamtools.btofmc.TestUtil.pairedMappedPrimarySamRecords;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

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
import com.hartwig.hmftools.bamtools.btofmc.TestUtil.ReadPairCollector;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class EvictingReadCacheTest
{
    private static class TestModule extends AbstractModule
    {
        public TestModule()
        {
        }

        @Override
        public void configure()
        {
            bind(EvictingReadCache.class).in(Singleton.class);

            bind(ReadPairCollector.class).in(Singleton.class);
            bind(new TypeLiteral<BiConsumer<SAMRecord, SAMRecord>>() {}).to(ReadPairCollector.class).in(Singleton.class);
        }
    }

    @Test
    public void testSizeOnEmpty()
    {
        Injector injector = Guice.createInjector(Stage.PRODUCTION, new TestModule());
        EvictingReadCache readCache = injector.getInstance(EvictingReadCache.class);
        ReadPairCollector readPairCollector = injector.getInstance(ReadPairCollector.class);

        assertEquals(0, readCache.size());
        assertTrue(readCache.isEmpty());
        assertEquals(0, readCache.addCount());
        assertEquals(0, readCache.evictionCount());
        assertTrue(readPairCollector.ReadPairs.isEmpty());
    }

    @Test
    public void testSizeOnAfterAddingSingleRead()
    {
        Injector injector = Guice.createInjector(Stage.PRODUCTION, new TestModule());
        EvictingReadCache readCache = injector.getInstance(EvictingReadCache.class);
        ReadPairCollector readPairCollector = injector.getInstance(ReadPairCollector.class);

        SAMRecord read = mappedPrimarySamRecord("READ_01", CHR_1, 1);
        readCache.add(read);

        assertEquals(1, readCache.size());
        assertFalse(readCache.isEmpty());
        assertEquals(1, readCache.addCount());
        assertEquals(0, readCache.evictionCount());
        assertTrue(readPairCollector.ReadPairs.isEmpty());
    }

    @Test
    public void testEvictionWithSingleBucket()
    {
        ReadPairCollector readPairCollector = new ReadPairCollector();
        EvictingReadCache readCache = new EvictingReadCache(readPairCollector, 1, false);

        readCache.add(mappedPrimarySamRecord("READ_01", CHR_1, 1));
        readCache.add(mappedPrimarySamRecord("READ_02", CHR_1, 1));

        assertEquals(2, readCache.size());
        assertFalse(readCache.isEmpty());
        assertEquals(2, readCache.addCount());
        assertEquals(1, readCache.evictionCount());
        assertTrue(readPairCollector.ReadPairs.isEmpty());
    }

    @Test
    public void testSizeOnAfterAddingPairedReadsWithoutEviction()
    {
        Injector injector = Guice.createInjector(Stage.PRODUCTION, new TestModule());
        EvictingReadCache readCache = injector.getInstance(EvictingReadCache.class);
        ReadPairCollector readPairCollector = injector.getInstance(ReadPairCollector.class);

        Pair<SAMRecord, SAMRecord> readPair = pairedMappedPrimarySamRecords("READ_01", CHR_1, 1, CHR_1, 1);
        readCache.add(readPair.getLeft());
        readCache.add(readPair.getRight());

        assertEquals(0, readCache.size());
        assertTrue(readCache.isEmpty());
        assertEquals(1, readCache.addCount());
        assertEquals(0, readCache.evictionCount());
        assertEquals(1, readPairCollector.ReadPairs.size());
    }

    @Test
    public void testAddingTwoReadPairsOutOfOrderWithConstantHashCode()
    {
        ReadPairCollector readPairCollector = new ReadPairCollector();
        EvictingReadCache readCache = new EvictingReadCache(readPairCollector, 2, false, record -> 0);

        Pair<SAMRecord, SAMRecord> readPair1 = pairedMappedPrimarySamRecords("READ_01", CHR_1, 1, CHR_1, 1);
        Pair<SAMRecord, SAMRecord> readPair2 = pairedMappedPrimarySamRecords("READ_02", CHR_1, 1, CHR_1, 1);

        readCache.add(readPair1.getLeft());
        readCache.add(readPair2.getLeft());
        readCache.add(readPair1.getRight());
        readCache.add(readPair2.getRight());

        assertEquals(4, readCache.size());
        assertFalse(readCache.isEmpty());
        assertEquals(4, readCache.addCount());
        assertEquals(3, readCache.evictionCount());
        assertTrue(readPairCollector.ReadPairs.isEmpty());

        readCache.flush();
        assertEquals(0, readCache.size());
        assertTrue(readCache.isEmpty());
        assertEquals(4, readCache.addCount());
        assertEquals(3, readCache.evictionCount());
        assertEquals(2, readPairCollector.ReadPairs.size());
    }

    @Test
    public void testConstantNegativeHashCode()
    {
        ReadPairCollector readPairCollector = new ReadPairCollector();
        EvictingReadCache readCache = new EvictingReadCache(readPairCollector, 32, false, record -> -1);

        SAMRecord read = mappedPrimarySamRecord("READ_01", CHR_1, 1);
        readCache.add(read);

        assertEquals(1, readCache.size());
        assertFalse(readCache.isEmpty());
        assertEquals(1, readCache.addCount());
        assertEquals(0, readCache.evictionCount());
        assertTrue(readPairCollector.ReadPairs.isEmpty());
    }

    @Test
    public void testForceResize()
    {
        ReadPairCollector readPairCollector = new ReadPairCollector();
        EvictingReadCache readCache = new EvictingReadCache(readPairCollector, 1);

        readCache.add(mappedPrimarySamRecord("READ_01", CHR_1, 1));
        readCache.add(mappedPrimarySamRecord("READ_02", CHR_1, 1));
        readCache.add(mappedPrimarySamRecord("READ_03", CHR_1, 1));
        readCache.add(mappedPrimarySamRecord("READ_04", CHR_1, 1));

        assertEquals(4, readCache.size());
        assertFalse(readCache.isEmpty());
        assertEquals(4, readCache.addCount());
        assertTrue(readPairCollector.ReadPairs.isEmpty());
    }

    @Test
    public void testOrphanReadWithLowestSortOrder()
    {
        ReadPairCollector readPairCollector = new ReadPairCollector();
        EvictingReadCache readCache = new EvictingReadCache(readPairCollector, 1, false);

        SAMRecord orphanedRead = mappedPrimarySamRecord("READ_01", CHR_1, 1);
        Pair<SAMRecord, SAMRecord> readPair = pairedMappedPrimarySamRecords("READ_02", CHR_1, 1, CHR_1, 1);

        readCache.add(readPair.getLeft());
        readCache.add(orphanedRead);
        readCache.add(readPair.getRight());

        readCache.flush();

        assertEquals(1, readCache.size());
        assertFalse(readCache.isEmpty());
        assertEquals(3, readCache.addCount());
        assertEquals(2, readCache.evictionCount());
        assertEquals(1, readPairCollector.ReadPairs.size());
    }

    @Test
    public void testOrphanReadWithHighestSortOrder()
    {
        ReadPairCollector readPairCollector = new ReadPairCollector();
        EvictingReadCache readCache = new EvictingReadCache(readPairCollector, 1, false);

        SAMRecord orphanedRead = mappedPrimarySamRecord("READ_02", CHR_1, 1);
        Pair<SAMRecord, SAMRecord> readPair = pairedMappedPrimarySamRecords("READ_01", CHR_1, 1, CHR_1, 1);

        readCache.add(readPair.getLeft());
        readCache.add(orphanedRead);
        readCache.add(readPair.getRight());

        readCache.flush();

        assertEquals(1, readCache.size());
        assertFalse(readCache.isEmpty());
        assertEquals(3, readCache.addCount());
        assertEquals(2, readCache.evictionCount());
        assertEquals(1, readPairCollector.ReadPairs.size());
    }
}
