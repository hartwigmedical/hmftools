package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.LiftBackScenario.ReadRole.MATE;
import static com.hartwig.hmftools.tars.liftback.LiftBackScenario.ReadRole.PRIMARY;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TX_CONTIG;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.nameGroupedHeader;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.pairedRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.supplementaryConfig;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.supplementaryRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.threeExonContig;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.writeBam;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.stream.Collectors;

import com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;
import com.hartwig.hmftools.tars.liftback.merge.SupplementaryMerger.RejectReason;
import com.hartwig.hmftools.tars.liftback.select.AlignmentSelector;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;

public class LiftBackPipelineTest
{
    @Rule
    public TemporaryFolder mFolder = new TemporaryFolder();

    @Test(timeout = 30_000)
    public void testBamOutputIsValidAndIndependentOfWorkerCount() throws Exception
    {
        String readBases = "ACGT".repeat(25);
        File ref = writeReference();

        SAMFileHeader inputHeader = nameGroupedHeader(
                new SAMSequenceRecord(CHR_1, 2_000), new SAMSequenceRecord(TX_CONTIG, 250));
        File input = mFolder.newFile("input.bam");
        List<SAMRecord> records = inputRecords(readBases);
        writeBam(input, inputHeader, records);

        List<SAMRecord> singleWorker = runPipeline(input, ref, 1, records.size());
        List<SAMRecord> multipleWorkers = runPipeline(input, ref, 3, records.size());
        assertEquals(sortedSamStrings(singleWorker), sortedSamStrings(multipleWorkers));

        LiftBackScenario.Result result = new LiftBackScenario.Result(multipleWorkers);
        result.assertEmittedCount(28);
        result.assertLifted("merge", PRIMARY, CHR_1, 150, "50M100N50M")
                .assertSuppCount("merge", 0)
                .assertTag("merge", PRIMARY, "AS", 100)
                .assertTag("merge", PRIMARY, "NM", 0)
                .assertTag("merge", PRIMARY, "SA", null)
                .assertTag("merge", MATE, "MC", "50M100N50M");
        result.assertUnmapped("lowScore", PRIMARY)
                .assertMapQuality("lowScore", PRIMARY, 0)
                .assertTag("lowScore", PRIMARY, "UM", CHR_1 + ":900")
                .assertLifted("lowScore", MATE, CHR_1, 800, "50M")
                .assertMateUnmapped("lowScore", MATE)
                .assertProperPair("lowScore", MATE, false)
                .assertTag("lowScore", MATE, "MC", null);
        assertEquals("*", result.record("lowScore", PRIMARY).getCigarString());
        for(int i = 0; i < 12; ++i)
        {
            String name = "fragment" + i;
            result.assertLifted(name, PRIMARY, CHR_1, 150, "50M100N50M")
                    .assertLifted(name, MATE, CHR_1, 800, "50M")
                    .assertMate(name, PRIMARY, CHR_1, 800)
                    .assertMate(name, MATE, CHR_1, 150)
                    .assertMapQuality(name, PRIMARY, 60)
                    .assertTag(name, PRIMARY, "NM", 0)
                    .assertTag(name, PRIMARY, "AS", 100)
                    .assertTag(name, PRIMARY, "MD", null)
                    .assertTag(name, PRIMARY, "MC", "50M")
                    .assertTag(name, MATE, "MC", "50M100N50M");
            SAMRecord first = result.record(name, PRIMARY);
            SAMRecord second = result.record(name, MATE);
            assertEquals(700, first.getInferredInsertSize());
            assertEquals(-700, second.getInferredInsertSize());
            assertEquals(readBases, first.getReadString());
            assertEquals("I".repeat(100), first.getBaseQualityString());
            assertFalse(first.getReadNegativeStrandFlag());
            assertTrue(first.getMateNegativeStrandFlag());
            assertTrue(second.getReadNegativeStrandFlag());
            assertFalse(second.getMateNegativeStrandFlag());
        }
    }

    @Test(timeout = 60_000)
    public void testLargeMixedInputPreservesRecordsAndCountersAcrossWorkersAndRepeatedRuns() throws Exception
    {
        File ref = writeReference();
        File input = writeMixedInput(1_000);
        PipelineOutput control = runPipeline(input, ref, 1, 1, 9_000);
        assertEquals(8_000, control.Records.size());
        assertEquals(Long.valueOf(1_000), control.Counters.get("SupplementariesAbsorbed"));
        assertEquals(Long.valueOf(1_000), control.Counters.get("UnmappedLowAlignmentScore"));
        Set<Integer> tieStarts = new HashSet<>();
        for(SAMRecord record : control.Records)
        {
            if(record.getReadName().endsWith("_tie") && record.getFirstOfPairFlag())
            {
                assertEquals("random ties keep MAPQ zero", 0, record.getMappingQuality());
                tieStarts.add(record.getAlignmentStart());
            }
        }
        assertEquals("read-seeded ties must exercise both placements", Set.of(1_000, 1_200), tieStarts);

        for(int workers : List.of(1, 2, 4))
        {
            for(int repeat = 0; repeat < 2; ++repeat)
            {
                PipelineOutput output = runPipeline(input, ref, workers, workers, 9_000);
                assertEquals("record contents for workers=" + workers + " repeat=" + repeat,
                        sortedSamStrings(control.Records), sortedSamStrings(output.Records));
                assertEquals("counters for workers=" + workers + " repeat=" + repeat, control.Counters, output.Counters);
            }
        }
    }

    File writeMixedInput(final int repetitions) throws Exception
    {
        List<SAMRecord> original = inputRecords("ACGT".repeat(25));
        List<SAMRecord> templates = new ArrayList<>(original.subList(0, 2));
        templates.addAll(original.subList(original.size() - 5, original.size()));
        SAMRecord first = pairedRecord("tie", true, CHR_1, 1_000, "50M", false);
        SAMRecord second = pairedRecord("tie", false, CHR_1, 1_050, "50M", true);
        first.setAttribute("XA", CHR_1 + ",+1200,50M,0;");
        second.setAttribute("XA", CHR_1 + ",-1250,50M,0;");
        for(SAMRecord record : List.of(first, second))
        {
            record.setMappingQuality(0);
            record.setReadString("A".repeat(50));
            record.setBaseQualityString("I".repeat(50));
            record.setAttribute("AS", 50);
            record.setAttribute("XS", 50);
        }
        templates.addAll(List.of(first, second));
        List<SAMRecord> records = new ArrayList<>();
        for(int i = 0; i < repetitions; ++i)
        {
            for(SAMRecord template : templates)
            {
                SAMRecord record = template.deepCopy();
                record.setReadName("repeat" + i + "_" + template.getReadName());
                records.add(record);
            }
        }
        File input = mFolder.newFile("mixed.bam");
        writeBam(input, nameGroupedHeader(new SAMSequenceRecord(CHR_1, 2_000), new SAMSequenceRecord(TX_CONTIG, 250)), records);
        return input;
    }

    File writeReference() throws Exception
    {
        String readBases = "ACGT".repeat(25);
        TestGenome genome = new TestGenome().with(CHR_1, 2_000, 'A')
                .set(CHR_1, 150, readBases.substring(0, 50))
                .set(CHR_1, 300, readBases.substring(50))
                .set(CHR_1, 800, "C".repeat(50));
        File ref = mFolder.newFile("genome.fa");
        Files.writeString(ref.toPath(), ">" + CHR_1 + "\n" + genome.asRefGenome().getBaseString(CHR_1, 1, 2_000) + "\n",
                StandardCharsets.US_ASCII);
        FastaSequenceIndexCreator.create(ref.toPath(), true);
        Files.writeString(mFolder.newFile("genome.dict").toPath(),
                "@HD\tVN:1.6\n@SQ\tSN:" + CHR_1 + "\tLN:2000\n", StandardCharsets.US_ASCII);

        return ref;
    }

    private static List<SAMRecord> inputRecords(final String readBases)
    {
        List<SAMRecord> records = new ArrayList<>();
        for(int i = 0; i < 12; ++i)
        {
            String name = "fragment" + i;
            SAMRecord first = pairedRecord(name, true, TX_CONTIG, 51, "100M", false);
            first.setReadString(readBases);
            first.setBaseQualityString("I".repeat(100));
            first.setMateReferenceName(CHR_1);
            first.setMateAlignmentStart(800);
            first.setMateNegativeStrandFlag(true);
            first.setAttribute("MD", "100");
            first.setAttribute("AS", 100);

            SAMRecord second = pairedRecord(name, false, CHR_1, 800, "50M", true);
            second.setReadString("C".repeat(50));
            second.setBaseQualityString("I".repeat(50));
            second.setMateReferenceName(TX_CONTIG);
            second.setMateAlignmentStart(51);
            second.setAttribute("AS", 50);
            records.add(first);
            records.add(second);
        }
        SAMRecord primary = pairedRecord("merge", true, CHR_1, 150, "50M50S", false);
        SAMRecord supplementary = supplementaryRecord("merge", CHR_1, 300, "50S50M");
        SAMRecord mate = pairedRecord("merge", false, CHR_1, 800, "50M", true);
        for(SAMRecord record : List.of(primary, supplementary))
        {
            record.setReadString(readBases);
            record.setBaseQualityString("I".repeat(100));
            record.setAttribute("AS", 50);
        }
        mate.setReadString("C".repeat(50));
        mate.setBaseQualityString("I".repeat(50));
        mate.setAttribute("AS", 50);
        records.addAll(List.of(primary, supplementary, mate));

        SAMRecord lowScore = pairedRecord("lowScore", true, CHR_1, 900, "10M", false);
        lowScore.setReadString("A".repeat(10));
        lowScore.setBaseQualityString("I".repeat(10));
        lowScore.setAttribute("AS", 10);
        SAMRecord lowScoreMate = mate.deepCopy();
        lowScoreMate.setReadName("lowScore");
        records.addAll(List.of(lowScore, lowScoreMate));
        return records;
    }

    private List<SAMRecord> runPipeline(final File input, final File reference, final int workerCount, final int expectedRecords)
            throws Exception
    {
        return runPipeline(input, reference, workerCount, 4, expectedRecords).Records;
    }

    PipelineOutput runPipeline(
            final File input, final File reference, final int workerCount, final int shardCount, final int expectedRecords)
            throws Exception
    {
        File outputDirectory = mFolder.newFolder();
        SAMFileHeader outputHeader = nameGroupedHeader(new SAMSequenceRecord(CHR_1, 2_000));
        LiftBackResources resources = new LiftBackResources(
                new AlignmentSelector(List.of(threeExonContig())),
                EnsemblAnnotationIndex.fromContigEntries(List.of(threeExonContig())),
                reference.getAbsolutePath(), supplementaryConfig(), null);
        BlockingQueue<List<SAMRecord>> queue = new LinkedBlockingQueue<>(2);
        List<LiftBackWorker> workers = new ArrayList<>();
        List<File> shards = new ArrayList<>();
        for(int i = 0; i < workerCount; ++i)
        {
            File shard = new File(outputDirectory, "output_" + i + ".bam");
            shards.add(shard);
            LiftBackWorker worker = new LiftBackWorker(queue, resources, outputHeader, shard.getAbsolutePath());
            worker.setDaemon(true);
            workers.add(worker);
            worker.start();
        }
        ShardedChunkProducer producer = new ShardedChunkProducer(
                List.of(input.getAbsolutePath()), reference.getAbsolutePath(), queue, workerCount, 128, shardCount);
        producer.setDaemon(true);
        producer.start();
        producer.join(15_000);
        assertFalse("producer did not terminate", producer.isAlive());

        long recordsSeen = 0;
        LiftBackStats totals = new LiftBackStats();
        for(LiftBackWorker worker : workers)
        {
            worker.join(15_000);
            assertFalse("worker did not terminate", worker.isAlive());
            recordsSeen += worker.stats().RecordsSeen;
            totals.add(worker.stats());
        }
        assertEquals(expectedRecords, recordsSeen);

        List<SAMRecord> emitted = new ArrayList<>();
        for(File shard : shards)
        {
            try(SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.STRICT).open(shard))
            {
                assertEquals(1, reader.getFileHeader().getSequenceDictionary().size());
                assertNull(reader.getFileHeader().getSequenceDictionary().getSequence(TX_CONTIG));
                for(SAMRecord record : reader)
                {
                    assertNull(record.getReadName() + " invalid SAM record", record.isValid());
                    emitted.add(record);
                }
            }
        }
        Map<String, Long> counters = new LinkedHashMap<>();
        for(Field field : LiftBackStats.class.getFields())
        {
            counters.put(field.getName(), field.getLong(totals));
        }
        for(RejectReason reason : RejectReason.values())
        {
            counters.put("Rejected_" + reason, totals.supplementaryMergeRejections(reason));
        }
        return new PipelineOutput(emitted, counters);
    }

    static class PipelineOutput
    {
        final List<SAMRecord> Records;
        final Map<String, Long> Counters;

        PipelineOutput(final List<SAMRecord> records, final Map<String, Long> counters)
        {
            Records = records;
            Counters = counters;
        }
    }

    private static List<String> sortedSamStrings(final List<SAMRecord> records)
    {
        return records.stream().map(SAMRecord::getSAMString).sorted().collect(Collectors.toList());
    }
}
