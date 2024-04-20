package com.hartwig.hmftools.bamtools.bamtofastq;

import static com.hartwig.hmftools.bamtools.bamtofastq.RegionTask.UNMAPPED_READS;
import static com.hartwig.hmftools.bamtools.bamtofastq.RegionTask.createChrRegion;
import static com.hartwig.hmftools.bamtools.bamtofastq.TestUtil.consensusSamRecord;
import static com.hartwig.hmftools.bamtools.bamtofastq.TestUtil.duplicateSamRecord;
import static com.hartwig.hmftools.bamtools.bamtofastq.TestUtil.mappedPrimarySamRecord;
import static com.hartwig.hmftools.bamtools.bamtofastq.TestUtil.supplementarySamRecord;
import static com.hartwig.hmftools.bamtools.bamtofastq.TestUtil.unmappedSamRecord;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;

import java.util.Collection;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Consumer;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;
import com.google.common.collect.Queues;
import com.google.inject.AbstractModule;
import com.google.inject.Guice;
import com.google.inject.Injector;
import com.google.inject.Provides;
import com.google.inject.Singleton;
import com.google.inject.Stage;
import com.google.inject.TypeLiteral;
import com.hartwig.hmftools.bamtools.bamtofastq.util.MockSamReader;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.commons.compress.utils.Lists;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class RegionTaskConsumerTest
{
    private static final ChrBaseRegion CHR_REGION = new ChrBaseRegion(CHR_1, 1_000, 2_000);
    private static final List<RegionTask> REGION_TASKS = List.of(UNMAPPED_READS, createChrRegion(CHR_REGION));

    private static class MockReadCache implements Consumer<SAMRecord>
    {
        public final List<SAMRecord> Reads;

        public MockReadCache()
        {
            Reads = Lists.newArrayList();
        }

        @Override
        public void accept(final SAMRecord samRecord)
        {
            Reads.add(samRecord);
        }
    }

    private static class TestModule extends AbstractModule
    {
        private final Collection<SAMRecord> mSamRecords;
        private final BamToFastqConfig mConfig;

        public TestModule(boolean keepConsensusReads, final Collection<SAMRecord> samRecords)
        {
            mSamRecords = samRecords;
            mConfig = new BamToFastqConfig(keepConsensusReads);
        }

        @Provides
        public RegionTaskConsumer provideRegionTaskConsumer(final RegionTaskConsumerFactory factory)
        {
            return factory.create();
        }

        @Provides
        @Singleton
        public Supplier<SamReader> provideSamReaderSupplier()
        {
            return () -> MockSamReader.create(mSamRecords);
        }

        @Override
        public void configure()
        {
            bind(BamToFastqConfig.class).toInstance(mConfig);
            bind(AtomicInteger.class).toInstance(new AtomicInteger());
            bind(new TypeLiteral<Queue<RegionTask>>() {}).toInstance(Queues.newArrayDeque(REGION_TASKS));

            bind(MockReadCache.class).in(Singleton.class);
            bind(new TypeLiteral<Consumer<SAMRecord>>() {}).to(MockReadCache.class).in(Singleton.class);

            bind(RegionTaskConsumerFactory.class).in(Singleton.class);
        }
    }

    private static class RegionTaskConsumerTestCase
    {
        private final boolean mKeepConsensusReads;
        private final Collection<SAMRecord> mReads;
        private final Multiset<String> mExpectedCachedReadNames;

        private RegionTaskConsumerTestCase(boolean keepConsensusReads, final Collection<SAMRecord> reads,
                final Multiset<String> expectedCachedReadNames)
        {
            mKeepConsensusReads = keepConsensusReads;
            mReads = reads;
            mExpectedCachedReadNames = expectedCachedReadNames;
        }

        public static void createAndCheck(boolean keepConsensusReads, final Collection<SAMRecord> reads,
                final Multiset<String> expectedCachedReadNames)
        {
            RegionTaskConsumerTestCase testCase = new RegionTaskConsumerTestCase(keepConsensusReads, reads, expectedCachedReadNames);
            testCase.check();
        }

        private void check()
        {
            Injector injector = Guice.createInjector(Stage.PRODUCTION, new TestModule(mKeepConsensusReads, mReads));
            RegionTaskConsumer consumer = injector.getInstance(RegionTaskConsumer.class);
            consumer.run();
            MockReadCache mockReadCache = injector.getInstance(MockReadCache.class);
            Multiset<String> actualCachedReadNames =
                    mockReadCache.Reads.stream().map(SAMRecord::getReadName).collect(Collectors.toCollection(() -> HashMultiset.create()));

            assertEquals(mExpectedCachedReadNames, actualCachedReadNames);
            assertEquals(mExpectedCachedReadNames.size(), consumer.readsProcessedCount());

            AtomicInteger totalRegionsProcessedCount = injector.getInstance(AtomicInteger.class);
            assertEquals(REGION_TASKS.size(), totalRegionsProcessedCount.get());
        }
    }

    @Test
    public void testUnmappedReadsAreProcessed()
    {
        SAMRecord read = unmappedSamRecord("READ_01");
        Multiset<String> expectedCachedReadNames = HashMultiset.create();
        expectedCachedReadNames.add("READ_01");

        RegionTaskConsumerTestCase.createAndCheck(false, List.of(read), expectedCachedReadNames);
    }

    @Test
    public void testMappedPrimaryReadsAreProcessed()
    {
        SAMRecord read = mappedPrimarySamRecord("READ_01", CHR_REGION.Chromosome, CHR_REGION.start());
        Multiset<String> expectedCachedReadNames = HashMultiset.create();
        expectedCachedReadNames.add("READ_01");

        RegionTaskConsumerTestCase.createAndCheck(false, List.of(read), expectedCachedReadNames);
    }

    @Test
    public void testDuplicateReadsAreProcessed()
    {
        SAMRecord read = duplicateSamRecord("READ_01", CHR_REGION.Chromosome, CHR_REGION.start());
        Multiset<String> expectedCachedReadNames = HashMultiset.create();
        expectedCachedReadNames.add("READ_01");

        RegionTaskConsumerTestCase.createAndCheck(false, List.of(read), expectedCachedReadNames);
    }

    @Test
    public void testSupplementaryReadsAreSkipped()
    {
        SAMRecord read = supplementarySamRecord("READ_01", CHR_REGION.Chromosome, CHR_REGION.start());
        Multiset<String> expectedCachedReadNames = HashMultiset.create();

        RegionTaskConsumerTestCase.createAndCheck(false, List.of(read), expectedCachedReadNames);
    }

    @Test
    public void testConsensusReads()
    {
        // Keep consensus reads.
        SAMRecord read = consensusSamRecord("READ_01", CHR_REGION.Chromosome, CHR_REGION.start());
        Multiset<String> expectedCachedReadNames = HashMultiset.create();
        expectedCachedReadNames.add("READ_01");

        RegionTaskConsumerTestCase.createAndCheck(true, List.of(read), expectedCachedReadNames);

        // Skip consensus reads.
        expectedCachedReadNames = HashMultiset.create();

        RegionTaskConsumerTestCase.createAndCheck(false, List.of(read), expectedCachedReadNames);
    }

    @Test
    public void testReadOverlappingButAlignmentStartNotInRegionIsSkipped()
    {
        SAMRecord read = mappedPrimarySamRecord("READ_01", CHR_REGION.Chromosome, CHR_REGION.start() - 1);
        Multiset<String> expectedCachedReadNames = HashMultiset.create();

        RegionTaskConsumerTestCase.createAndCheck(false, List.of(read), expectedCachedReadNames);
    }

    @Test
    public void testMergeStats()
    {
        RegionTaskConsumer consumer1 = new RegionTaskConsumer(2);
        RegionTaskConsumer consumer2 = new RegionTaskConsumer(3);
        consumer1.mergeStats(consumer2);

        assertEquals(5, consumer1.readsProcessedCount());
    }
}
