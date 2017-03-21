package com.hartwig.hmftools.fastqstats;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.net.URL;

import com.google.common.io.Resources;

import org.junit.Test;

public class FastqStatsTest {
    // 100% of baseCalls have a quality score >30
    @Test
    public void computesCorrectStatsFor100PercentDir() throws IOException, InterruptedException {
        final URL path = Resources.getResource("DeterminedBaseCalls-100");
        final File dir = new File(path.getPath());
        final FastqTracker tracker = FastqStats.processDir(dir);
        assertEquals(300, tracker.flowcell().yield());
        final double q30Percentage = tracker.flowcell().q30() * 100.0 / tracker.flowcell().yield();
        assertEquals(100.0, q30Percentage, 0.00001);
        assertEquals(300, tracker.lane("L001").yield());
        assertEquals(3, tracker.samples().size());
        final double laneQ30Percentage = tracker.lane("L001").q30() * 100.0 / tracker.lane("L001").yield();
        assertEquals(100.0, laneQ30Percentage, 0.00001);
    }

    // Sample Z has 50% Q30, others 100%
    @Test
    public void computesCorrectStatsForIndividualSamples() throws IOException, InterruptedException {
        final URL path = Resources.getResource("DeterminedBaseCalls-100&50");
        final File dir = new File(path.getPath());
        final FastqTracker tracker = FastqStats.processDir(dir);
        assertEquals(300, tracker.flowcell().yield());
        final double q30Percentage = tracker.flowcell().q30() * 100.0 / tracker.flowcell().yield();
        assertEquals(83.33, q30Percentage, 0.01);
        assertEquals(300, tracker.lane("L001").yield());
        final double laneQ30Percentage = tracker.lane("L001").q30() * 100.0 / tracker.lane("L001").yield();
        assertEquals(83.33, laneQ30Percentage, 0.01);
        assertEquals(3, tracker.samples().size());
        final double sampleXQ30 = tracker.sample("TESTX").q30() * 100.0 / tracker.sample("TESTX").yield();
        final double sampleYQ30 = tracker.sample("TESTY").q30() * 100.0 / tracker.sample("TESTY").yield();
        final double sampleZQ30 = tracker.sample("TESTZ").q30() * 100.0 / tracker.sample("TESTZ").yield();
        assertEquals(100.0, sampleXQ30, 0.00001);
        assertEquals(100.0, sampleYQ30, 0.00001);
        assertEquals(50.0, sampleZQ30, 0.00001);
    }

    // 1 Sample (X), 3 Lanes; Lane 3 has 50% Q30, others 100%
    @Test
    public void computesCorrectStatsForIndividualLanes() throws IOException, InterruptedException {
        final URL path = Resources.getResource("DeterminedBaseCallsMultiLane-100&50");
        final File dir = new File(path.getPath());
        final FastqTracker tracker = FastqStats.processDir(dir);
        assertEquals(300, tracker.flowcell().yield());
        final double q30Percentage = tracker.flowcell().q30() * 100.0 / tracker.flowcell().yield();
        assertEquals(83.33, q30Percentage, 0.01);
        assertEquals(1, tracker.samples().size());
        assertEquals(3, tracker.lanes().size());
        assertEquals(100, tracker.lane("L001").yield());
        assertEquals(100, tracker.lane("L002").yield());
        assertEquals(100, tracker.lane("L003").yield());
        assertEquals(100, tracker.lane("L001").q30() * 100.0 / tracker.lane("L001").yield(), 0.01);
        assertEquals(100, tracker.lane("L002").q30() * 100.0 / tracker.lane("L002").yield(), 0.01);
        assertEquals(50, tracker.lane("L003").q30() * 100.0 / tracker.lane("L003").yield(), 0.01);
    }

    // 50% of baseCalls belong to undetermined sample
    @Test
    public void computesCorrectPercentageOfUndetermined50() throws IOException, InterruptedException {
        final URL path = Resources.getResource("50UndeterminedBaseCalls");
        final File dir = new File(path.getPath());
        final FastqTracker tracker = FastqStats.processDir(dir);
        assertEquals(600, tracker.flowcell().yield());
        final double q30Percentage = tracker.flowcell().q30() * 100.0 / tracker.flowcell().yield();
        assertEquals(100.0, q30Percentage, 0.00001);
        assertEquals(3, tracker.samples().size());
        final double undeterminedPercentage = tracker.undetermined().yield() * 100.0 / tracker.flowcell().yield();
        assertEquals(50.0, undeterminedPercentage, 0.00001);
        assertEquals(600, tracker.lane("L001").yield());
    }

    // 100% of baseCalls belong to undetermined sample
    @Test
    public void computesCorrectPercentageOfUndetermined100() throws IOException, InterruptedException {
        final URL path = Resources.getResource("100UndeterminedBaseCalls");
        final File dir = new File(path.getPath());
        final FastqTracker tracker = FastqStats.processDir(dir);
        assertEquals(300, tracker.flowcell().yield());
        assertEquals(100.0, tracker.flowcell().q30() * 100.0 / tracker.flowcell().yield(), 0.00001);
        assertEquals(0, tracker.samples().size());
        final double undeterminedPercentage = tracker.undetermined().yield() * 100.0 / tracker.flowcell().yield();
        assertEquals(100.0, undeterminedPercentage, 0.00001);
        assertEquals(300, tracker.lane("L002").yield());
    }

    // 50% of baseCalls are undetermined and all undetermined have Q30 of 0
    @Test
    public void computesCorrectStatsOfUndeterminedLane() throws IOException, InterruptedException {
        final URL path = Resources.getResource("UndeterminedBaseCalls0Q30");
        final File dir = new File(path.getPath());
        final FastqTracker tracker = FastqStats.processDir(dir);
        assertEquals(600, tracker.flowcell().yield());
        assertEquals(50.0, tracker.flowcell().q30() * 100.0 / tracker.flowcell().yield(), 0.00001);
        final double undeterminedPercentage = tracker.undetermined().yield() * 100.0 / tracker.flowcell().yield();
        assertEquals(50.0, undeterminedPercentage, 0.00001);
        assertEquals(3, tracker.samples().size());
        final double undeterminedQ30 = tracker.undetermined().q30() * 100.0 / tracker.undetermined().yield();
        assertEquals(0.0, undeterminedQ30, 0.00001);
        assertEquals(300, tracker.lane("L002").yield());
    }

    // 100% of baseCalls have a quality score >30, but sample X is corrupt (has .fastq.gz extension, but is plain text).
    @Test
    public void skipsFailedSampleFile() throws IOException, InterruptedException {
        final URL path = Resources.getResource("CorruptSampleXFile");
        final File dir = new File(path.getPath());
        final FastqTracker tracker = FastqStats.processDir(dir);
        assertEquals(200, tracker.flowcell().yield());
        final double q30Percentage = tracker.flowcell().q30() * 100.0 / tracker.flowcell().yield();
        assertEquals(100.0, q30Percentage, 0.00001);
        assertEquals(200, tracker.lane("L001").yield());
        assertEquals(2, tracker.samples().size());
        final double laneQ30Percentage = tracker.lane("L001").q30() * 100.0 / tracker.lane("L001").yield();
        assertEquals(100.0, laneQ30Percentage, 0.00001);
    }

    // 100% of baseCalls belong to undetermined sample, but 1 file is corrupt (has .fastq.gz extension, but is plain text).
    @Test
    public void skipsFailedUndeterminedSampleFile() throws IOException, InterruptedException {
        final URL path = Resources.getResource("CorruptUndeterminedSampleFile");
        final File dir = new File(path.getPath());
        final FastqTracker tracker = FastqStats.processDir(dir);
        assertEquals(200, tracker.flowcell().yield());
        assertEquals(100.0, tracker.flowcell().q30() * 100.0 / tracker.flowcell().yield(), 0.00001);
        assertEquals(0, tracker.samples().size());
        final double undeterminedPercentage = tracker.undetermined().yield() * 100.0 / tracker.flowcell().yield();
        assertEquals(100.0, undeterminedPercentage, 0.00001);
        assertEquals(200, tracker.lane("L002").yield());
    }

    @Test
    public void computesCorrectStatsOfSingleZippedFile() throws IOException {
        final URL path = Resources.getResource("fastq" + File.separator + "q30-10_Flowcell_S1_L001_R1_001.fastq.gz");
        final FastqTracker tracker = FastqStats.processFile(path.getPath());
        assertEquals(100, tracker.flowcell().yield());
        assertEquals(10.0, tracker.flowcell().q30(), 0.00001);
    }

    @Test(expected = IOException.class)
    public void processFileThrowsForUnrecognizedFormat() throws IOException {
        final URL path = Resources.getResource("fastq" + File.separator + "q30-10.err");
        FastqStats.processFile(path.getPath());
    }

}
