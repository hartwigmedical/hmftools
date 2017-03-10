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
    public void computesCorrectStatsFor100PercentDir() throws IOException {
        final URL path = Resources.getResource("DeterminedBaseCalls-100");
        final File dir = new File(path.getPath());
        final FastqTracker tracker = FastqStats.processDir(dir);
        assertEquals(300, tracker.getFlowcellData().getYield());
        final double q30Percentage = tracker.getFlowcellData().getQ30() * 100.0 / tracker.getFlowcellData().getYield();
        assertEquals(100.0, q30Percentage, 0.00001);
        assertEquals(300, tracker.getLaneData("L001").getYield());
        final double laneQ30Percentage =
                tracker.getLaneData("L001").getQ30() * 100.0 / tracker.getLaneData("L001").getYield();
        assertEquals(100.0, laneQ30Percentage, 0.00001);
    }

    // Sample Z has 50% Q30, others 100%
    @Test
    public void computesCorrectStatsForIndividualSamples() throws IOException {
        final URL path = Resources.getResource("DeterminedBaseCalls-100&50");
        final File dir = new File(path.getPath());
        final FastqTracker tracker = FastqStats.processDir(dir);
        assertEquals(300, tracker.getFlowcellData().getYield());
        final double q30Percentage = tracker.getFlowcellData().getQ30() * 100.0 / tracker.getFlowcellData().getYield();
        assertEquals(83.33, q30Percentage, 0.01);
        assertEquals(300, tracker.getLaneData("L001").getYield());
        final double laneQ30Percentage =
                tracker.getLaneData("L001").getQ30() * 100.0 / tracker.getLaneData("L001").getYield();
        assertEquals(83.33, laneQ30Percentage, 0.01);
        final double sampleXQ30 =
                tracker.getSampleData("TESTX").getQ30() * 100.0 / tracker.getSampleData("TESTX").getYield();
        final double sampleYQ30 =
                tracker.getSampleData("TESTY").getQ30() * 100.0 / tracker.getSampleData("TESTY").getYield();
        final double sampleZQ30 =
                tracker.getSampleData("TESTZ").getQ30() * 100.0 / tracker.getSampleData("TESTZ").getYield();
        assertEquals(100.0, sampleXQ30, 0.00001);
        assertEquals(100.0, sampleYQ30, 0.00001);
        assertEquals(50.0, sampleZQ30, 0.00001);
    }

    // 50% of baseCalls belong to undetermined sample
    @Test
    public void computesCorrectPercentageOfUndetermined50() throws IOException {
        final URL path = Resources.getResource("50UndeterminedBaseCalls");
        final File dir = new File(path.getPath());
        final FastqTracker tracker = FastqStats.processDir(dir);
        assertEquals(600, tracker.getFlowcellData().getYield());
        final double q30Percentage = tracker.getFlowcellData().getQ30() * 100.0 / tracker.getFlowcellData().getYield();
        assertEquals(100.0, q30Percentage, 0.00001);
        final double undeterminedPercentage =
                tracker.getUndeterminedData().getYield() * 100.0 / tracker.getFlowcellData().getYield();
        assertEquals(50.0, undeterminedPercentage, 0.00001);
        assertEquals(600, tracker.getLaneData("L001").getYield());
    }

    // 100% of baseCalls belong to undetermined sample
    @Test
    public void computesCorrectPercentageOfUndetermined100() throws IOException {
        final URL path = Resources.getResource("100UndeterminedBaseCalls");
        final File dir = new File(path.getPath());
        final FastqTracker tracker = FastqStats.processDir(dir);
        assertEquals(300, tracker.getFlowcellData().getYield());
        assertEquals(100.0, tracker.getFlowcellData().getQ30() * 100.0 / tracker.getFlowcellData().getYield(),
                0.00001);
        final double undeterminedPercentage =
                tracker.getUndeterminedData().getYield() * 100.0 / tracker.getFlowcellData().getYield();
        assertEquals(100.0, undeterminedPercentage, 0.00001);
        assertEquals(300, tracker.getLaneData("L002").getYield());
    }

    // 50% of baseCalls are undetermined and all undetermined have Q30 of 0
    @Test
    public void computesCorrectStatsOfUndeterminedLane() throws IOException {
        final URL path = Resources.getResource("UndeterminedBaseCalls0Q30");
        final File dir = new File(path.getPath());
        final FastqTracker tracker = FastqStats.processDir(dir);
        assertEquals(600, tracker.getFlowcellData().getYield());
        assertEquals(50.0, tracker.getFlowcellData().getQ30() * 100.0 / tracker.getFlowcellData().getYield(), 0.00001);
        final double undeterminedPercentage =
                tracker.getUndeterminedData().getYield() * 100.0 / tracker.getFlowcellData().getYield();
        assertEquals(50.0, undeterminedPercentage, 0.00001);
        final double undeterminedQ30 =
                tracker.getUndeterminedData().getQ30() * 100.0 / tracker.getUndeterminedData().getYield();
        assertEquals(0.0, undeterminedQ30, 0.00001);
        assertEquals(300, tracker.getLaneData("L002").getYield());
    }
}
