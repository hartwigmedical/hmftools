package com.hartwig.hmftools.fastqstats.test;

import static org.junit.Assert.assertEquals;
import java.io.File;
import java.io.IOException;
import com.hartwig.hmftools.fastqstats.FastqStats;
import com.hartwig.hmftools.fastqstats.FastqStatsRunner;
import com.hartwig.hmftools.fastqstats.FastqTracker;
import org.junit.Test;

public class FastqStatsTest {
    // 100% of baseCalls have a quality score >30
    @Test
    public void computesCorrectStatsFor100PercentDir() throws IOException{
        File dir = new File("./src/test/resources/DeterminedBaseCalls-100");
        FastqTracker tracker = FastqStats.processDir(dir);
        assertEquals(300, tracker.getFlowcellData().getYield());
        double q30Percentage = tracker.getFlowcellData().getQ30() * 100.0 / tracker.getFlowcellData().getYield();
        assertEquals(100.0, q30Percentage, 0.00001);
        assertEquals(300, tracker.getLaneData("L001").getYield());
        double laneQ30Percentage = tracker.getLaneData("L001").getQ30() * 100.0 / tracker.getLaneData("L001").getYield();
        assertEquals(100.0, laneQ30Percentage, 0.00001);
    }

    // Sample Z has 50% Q30, others 100%
    @Test
    public void computesCorrectStatsForIndividualSamples() throws IOException{
        File dir = new File("./src/test/resources/DeterminedBaseCalls-100&50");
        FastqTracker tracker = FastqStats.processDir(dir);
        assertEquals(300, tracker.getFlowcellData().getYield());
        double q30Percentage = tracker.getFlowcellData().getQ30() * 100.0 / tracker.getFlowcellData().getYield();
        assertEquals(83.33, q30Percentage, 0.01);
        assertEquals(300, tracker.getLaneData("L001").getYield());
        double laneQ30Percentage = tracker.getLaneData("L001").getQ30() * 100.0 / tracker.getLaneData("L001").getYield();
        assertEquals(83.33, laneQ30Percentage, 0.01);
        double sampleXQ30 = tracker.getSampleData("TESTX").getQ30() * 100.0 / tracker.getSampleData("TESTX").getYield();
        double sampleYQ30 = tracker.getSampleData("TESTY").getQ30() * 100.0 / tracker.getSampleData("TESTY").getYield();
        double sampleZQ30 = tracker.getSampleData("TESTZ").getQ30() * 100.0 / tracker.getSampleData("TESTZ").getYield();
        assertEquals(100.0, sampleXQ30, 0.00001);
        assertEquals(100.0, sampleYQ30, 0.00001);
        assertEquals(50.0, sampleZQ30, 0.00001);
    }

    // 50% of baseCalls belong to undetermined sample
    @Test
    public void computesCorrectPercentageOfUndetermined50() throws IOException{
        File dir = new File("./src/test/resources/50UndeterminedBaseCalls");
        FastqTracker tracker = FastqStats.processDir(dir);
        assertEquals(600, tracker.getFlowcellData().getYield());
        double q30Percentage = tracker.getFlowcellData().getQ30() * 100.0 / tracker.getFlowcellData().getYield();
        assertEquals(100.0, q30Percentage, 0.00001);
        double undeterminedPercentage = tracker.getUndeterminedData().getYield() * 100.0 / tracker.getFlowcellData().getYield();
        assertEquals(50.0, undeterminedPercentage, 0.00001);
        assertEquals(600, tracker.getLaneData("L001").getYield());
        FastqStatsRunner.printOutput(tracker);
    }

    // 100% of baseCalls belong to undetermined sample
    @Test
    public void computesCorrectPercentageOfUndetermined100() throws IOException{
        File dir = new File("./src/test/resources/100UndeterminedBaseCalls");
        FastqTracker tracker = FastqStats.processDir(dir);
        assertEquals(300, tracker.getFlowcellData().getYield());
        assertEquals(100.0, tracker.getFlowcellData().getQ30() * 100.0 / tracker.getFlowcellData().getYield(), 0.00001);
        double undeterminedPercentage = tracker.getUndeterminedData().getYield() * 100.0 / tracker.getFlowcellData().getYield();
        assertEquals(100.0, undeterminedPercentage, 0.00001);
        assertEquals(300, tracker.getLaneData("L002").getYield());
        FastqStatsRunner.printOutput(tracker);
    }

    // 50% of baseCalls are undetermined and all undetermined have Q30 of 0
    @Test
    public void computesCorrectStatsOfUndeterminedLane() throws IOException{
        File dir = new File("./src/test/resources/UndeterminedBaseCalls0Q30");
        FastqTracker tracker = FastqStats.processDir(dir);
        assertEquals(600, tracker.getFlowcellData().getYield());
        assertEquals(50.0, tracker.getFlowcellData().getQ30() * 100.0 / tracker.getFlowcellData().getYield(), 0.00001);
        double undeterminedPercentage = tracker.getUndeterminedData().getYield() * 100.0 / tracker.getFlowcellData().getYield();
        assertEquals(50.0, undeterminedPercentage, 0.00001);
        double undeterminedQ30 = tracker.getUndeterminedData().getQ30() * 100.0 / tracker.getUndeterminedData().getYield();
        assertEquals(0.0, undeterminedQ30, 0.00001);
        assertEquals(300, tracker.getLaneData("L002").getYield());
        FastqStatsRunner.printOutput(tracker);
    }


}
