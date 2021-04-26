package com.hartwig.hmftools.fastqstats;

import static com.hartwig.hmftools.fastqstats.FastqStats.getFastqsFromBaseCallsDir;
import static com.hartwig.hmftools.fastqstats.FastqStats.getFastqsFromDir;
import static com.hartwig.hmftools.fastqstats.FastqStats.getSingleFastq;
import static com.hartwig.hmftools.fastqstats.FastqStatsRunner.getDir;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.net.URL;

import com.google.common.collect.Multimap;
import com.google.common.io.Resources;

import org.junit.Test;

public class FastqStatsTest {

    // 100% of baseCalls have a quality score >30
    @Test
    public void computesCorrectStatsFor100PercentDir() throws IOException, InterruptedException {
        final URL path = Resources.getResource("DeterminedBaseCalls-100");
        final File baseCallsDir = new File(path.getPath());
        final Multimap<String, File> fastqsPerSample = getFastqsFromBaseCallsDir(baseCallsDir);
        final FastqTracker tracker = FastqStats.processFastqs(fastqsPerSample, FastqStatsRunner.getThreadCount(null));
        assertEquals(300, tracker.flowcell().yield());
        assertEquals(100.0, tracker.flowcell().q30Percentage(), 0.00001);
        assertEquals(300, tracker.lane("L001").yield());
        assertEquals(3, tracker.samples().size());
        assertEquals(100.0, tracker.lane("L001").q30Percentage(), 0.00001);
    }

    // Sample Z has 50% Q30, others 100%
    @Test
    public void computesCorrectStatsForIndividualSamples() throws IOException, InterruptedException {
        final URL path = Resources.getResource("DeterminedBaseCalls-100&50");
        final File baseCallsDir = new File(path.getPath());
        final Multimap<String, File> fastqsPerSample = getFastqsFromBaseCallsDir(baseCallsDir);
        final FastqTracker tracker = FastqStats.processFastqs(fastqsPerSample, FastqStatsRunner.getThreadCount(null));
        assertEquals(300, tracker.flowcell().yield());
        assertEquals(83.33, tracker.flowcell().q30Percentage(), 0.01);
        assertEquals(300, tracker.lane("L001").yield());
        assertEquals(83.33, tracker.lane("L001").q30Percentage(), 0.01);
        assertEquals(3, tracker.samples().size());
        assertEquals(100.0, tracker.sample("TESTX").q30Percentage(), 0.00001);
        assertEquals(100.0, tracker.sample("TESTY").q30Percentage(), 0.00001);
        assertEquals(50.0, tracker.sample("TESTZ").q30Percentage(), 0.00001);
        assertEquals(100, tracker.samples().get("TESTZ").get("L001").yield());
        assertEquals(50, tracker.samples().get("TESTZ").get("L001").q30());
    }

    // 1 Sample (X), 3 Lanes; Lane 3 has 50% Q30, others 100%
    @Test
    public void computesCorrectStatsForIndividualLanes() throws IOException, InterruptedException {
        final URL path = Resources.getResource("DeterminedBaseCallsMultiLane-100&50");
        final File baseCallsDir = new File(path.getPath());
        final Multimap<String, File> fastqsPerSample = getFastqsFromBaseCallsDir(baseCallsDir);
        final FastqTracker tracker = FastqStats.processFastqs(fastqsPerSample, FastqStatsRunner.getThreadCount(null));
        assertEquals(300, tracker.flowcell().yield());
        assertEquals(83.33, tracker.flowcell().q30Percentage(), 0.01);
        assertEquals(1, tracker.samples().size());
        assertEquals(3, tracker.lanes().size());
        assertEquals(100, tracker.lane("L001").yield());
        assertEquals(100, tracker.lane("L002").yield());
        assertEquals(100, tracker.lane("L003").yield());
        assertEquals(100, tracker.samples().get("TESTX").get("L003").yield());
        assertEquals(50, tracker.samples().get("TESTX").get("L003").q30());
        assertEquals(100, tracker.lane("L001").q30Percentage(), 0.01);
        assertEquals(100, tracker.lane("L002").q30Percentage(), 0.01);
        assertEquals(50, tracker.lane("L003").q30Percentage(), 0.01);
    }

    // 50% of baseCalls belong to undetermined sample
    @Test
    public void computesCorrectPercentageOfUndetermined50() throws IOException, InterruptedException {
        final URL path = Resources.getResource("50UndeterminedBaseCalls");
        final File baseCallsDir = new File(path.getPath());
        final Multimap<String, File> fastqsPerSample = getFastqsFromBaseCallsDir(baseCallsDir);
        final FastqTracker tracker = FastqStats.processFastqs(fastqsPerSample, FastqStatsRunner.getThreadCount(null));
        assertEquals(600, tracker.flowcell().yield());
        assertEquals(100.0, tracker.flowcell().q30Percentage(), 0.00001);
        assertEquals(4, tracker.samples().size());
        final double undeterminedPercentage = tracker.sample("Undetermined").yield() * 100.0 / tracker.flowcell().yield();
        assertEquals(50.0, undeterminedPercentage, 0.00001);
        assertEquals(600, tracker.lane("L001").yield());
        assertEquals(300, tracker.samples().get("Undetermined").get("L001").yield());
        assertEquals(300, tracker.samples().get("Undetermined").get("L001").q30());
    }

    // 100% of baseCalls belong to undetermined sample
    @Test
    public void computesCorrectPercentageOfUndetermined100() throws IOException, InterruptedException {
        final URL path = Resources.getResource("100UndeterminedBaseCalls");
        final File baseCallsDir = new File(path.getPath());
        final Multimap<String, File> fastqsPerSample = getFastqsFromBaseCallsDir(baseCallsDir);
        final FastqTracker tracker = FastqStats.processFastqs(fastqsPerSample, FastqStatsRunner.getThreadCount(null));
        assertEquals(300, tracker.flowcell().yield());
        assertEquals(100.0, tracker.flowcell().q30Percentage(), 0.00001);
        assertEquals(1, tracker.samples().size());
        final double undeterminedPercentage = tracker.sample("Undetermined").yield() * 100.0 / tracker.flowcell().yield();
        assertEquals(100.0, undeterminedPercentage, 0.00001);
        assertEquals(300, tracker.lane("L002").yield());
        assertEquals(300, tracker.samples().get("Undetermined").get("L002").yield());
        assertEquals(300, tracker.samples().get("Undetermined").get("L002").q30());
    }

    // 50% of baseCalls are undetermined and all undetermined have Q30 of 0
    @Test
    public void computesCorrectStatsOfUndeterminedLane() throws IOException, InterruptedException {
        final URL path = Resources.getResource("UndeterminedBaseCalls0Q30");
        final File baseCallsDir = new File(path.getPath());
        final Multimap<String, File> fastqsPerSample = getFastqsFromBaseCallsDir(baseCallsDir);
        final FastqTracker tracker = FastqStats.processFastqs(fastqsPerSample, FastqStatsRunner.getThreadCount(null));
        assertEquals(600, tracker.flowcell().yield());
        assertEquals(50.0, tracker.flowcell().q30Percentage(), 0.00001);
        final double undeterminedPercentage = tracker.sample("Undetermined").yield() * 100.0 / tracker.flowcell().yield();
        assertEquals(50.0, undeterminedPercentage, 0.00001);
        assertEquals(4, tracker.samples().size());
        assertEquals(0.0, tracker.sample("Undetermined").q30Percentage(), 0.00001);
        assertEquals(300, tracker.lane("L002").yield());
        assertEquals(300, tracker.samples().get("Undetermined").get("L002").yield());
        assertEquals(0, tracker.samples().get("Undetermined").get("L002").q30());
    }

    // 100% of baseCalls have a quality score >30, but sample X is corrupt (has .fastq.gz extension, but is plain text).
    @Test
    public void skipsFailedSampleFile() throws IOException, InterruptedException {
        final URL path = Resources.getResource("CorruptSampleXFile");
        final File baseCallsDir = new File(path.getPath());
        final Multimap<String, File> fastqsPerSample = getFastqsFromBaseCallsDir(baseCallsDir);
        final FastqTracker tracker = FastqStats.processFastqs(fastqsPerSample, FastqStatsRunner.getThreadCount(null));
        assertEquals(200, tracker.flowcell().yield());
        assertEquals(100.0, tracker.flowcell().q30Percentage(), 0.00001);
        assertEquals(200, tracker.lane("L001").yield());
        assertEquals(2, tracker.samples().size());
        final double laneQ30Percentage = tracker.lane("L001").q30Percentage();
        assertEquals(100.0, laneQ30Percentage, 0.00001);
    }

    // 100% of baseCalls belong to undetermined sample, but 1 file is corrupt (has .fastq.gz extension, but is plain text).
    @Test
    public void skipsFailedUndeterminedSampleFile() throws IOException, InterruptedException {
        final URL path = Resources.getResource("CorruptUndeterminedSampleFile");
        final File baseCallsDir = new File(path.getPath());
        final Multimap<String, File> fastqsPerSample = getFastqsFromBaseCallsDir(baseCallsDir);
        final FastqTracker tracker = FastqStats.processFastqs(fastqsPerSample, FastqStatsRunner.getThreadCount(null));
        assertEquals(200, tracker.flowcell().yield());
        assertEquals(100.0, tracker.flowcell().q30Percentage(), 0.00001);
        assertEquals(1, tracker.samples().size());
        final double undeterminedPercentage = tracker.sample("Undetermined").yield() * 100.0 / tracker.flowcell().yield();
        assertEquals(100.0, undeterminedPercentage, 0.00001);
        assertEquals(200, tracker.lane("L002").yield());
    }

    @Test
    public void computesCorrectStatsOfSingleZippedFile() throws InterruptedException {
        final URL path = Resources.getResource("fastq/q30-10_Flowcell_S1_L001_R1_001.fastq.gz");
        final Multimap<String, File> fastqsPerSample = getSingleFastq(path.getPath());
        final FastqTracker tracker = FastqStats.processFastqs(fastqsPerSample, 1);
        assertEquals(100, tracker.flowcell().yield());
        assertEquals(10.0, tracker.flowcell().q30(), 0.00001);
    }

    @Test
    public void skipsUnrecognizedFormat() throws InterruptedException {
        final URL path = Resources.getResource("fastq/q30-10.err");
        final Multimap<String, File> fastqsPerSample = getSingleFastq(path.getPath());
        final FastqTracker tracker = FastqStats.processFastqs(fastqsPerSample, 1);
        assertEquals(0, tracker.flowcell().yield());
        assertEquals(0, tracker.flowcell().q30(), 0.00001);
    }

    @Test
    public void computesCorrectStatsForFastqDir() throws IOException, InterruptedException {
        final URL path = Resources.getResource("100UndeterminedBaseCalls");
        final File fastqDir = getDir(path.getPath());
        final Multimap<String, File> fastqs = getFastqsFromDir(fastqDir);
        final FastqTracker tracker = FastqStats.processFastqs(fastqs, FastqStatsRunner.getThreadCount(null));
        assertEquals(300, tracker.flowcell().yield());
        assertEquals(100.0, tracker.flowcell().q30Percentage(), 0.00001);
        assertEquals(3, tracker.samples().size());
        assertEquals(300, tracker.lane("L002").yield());
        assertEquals(100, tracker.sample("UndeterminedF1_H7YRLADXX_S1_L002_R1_001.fastq").yield());
        assertEquals(100, tracker.sample("UndeterminedF2_H7YRLADXX_S1_L002_R1_001.fastq").yield());
        assertEquals(100, tracker.sample("UndeterminedF3_H7YRLADXX_S1_L002_R1_001.fastq").yield());
        assertEquals(100, tracker.sample("UndeterminedF1_H7YRLADXX_S1_L002_R1_001.fastq").q30());
        assertEquals(100, tracker.sample("UndeterminedF2_H7YRLADXX_S1_L002_R1_001.fastq").q30());
        assertEquals(100, tracker.sample("UndeterminedF3_H7YRLADXX_S1_L002_R1_001.fastq").q30());
    }
}
