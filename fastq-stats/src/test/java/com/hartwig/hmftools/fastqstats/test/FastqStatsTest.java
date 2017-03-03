package com.hartwig.hmftools.fastqstats.test;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.fastqstats.FastqStats;
import com.hartwig.hmftools.fastqstats.FastqTracker;
import com.hartwig.hmftools.fastqstats.TrackerKey;
import static com.hartwig.hmftools.fastqstats.TrackerType.*;
import org.junit.Test;

public class FastqStatsTest {
    @Test
    public void computesCorrectStatsFor100PercentDir() throws IOException{
        File dir = new File("./src/test/resources/DeterminedBaseCalls-100");
        FastqTracker tracker = FastqStats.processDir(dir);
        TrackerKey yield = new TrackerKey(Flowcell, "Yield");
        TrackerKey q30 = new TrackerKey(Flowcell, "Q30");
        assertEquals(300, tracker.get(yield));
        assertEquals(100.0, tracker.getPercentage(q30, yield), 0.00001);
        TrackerKey laneYield = new TrackerKey(Lane, "L001-Yield");
        TrackerKey laneQ30 = new TrackerKey(Lane, "L001-Q30");
        assertEquals(300, tracker.get(laneYield));
        assertEquals(100.0, tracker.getPercentage(laneQ30, laneYield), 0.00001);
    }

    @Test
    public void computesCorrectStatsOfUndeterminedDir() throws IOException{
        File dir = new File("./src/test/resources/50UndeterminedBaseCalls");
        FastqTracker tracker = FastqStats.processDir(dir);
        TrackerKey yield = new TrackerKey(Flowcell, "Yield");
        TrackerKey q30 = new TrackerKey(Flowcell, "Q30");
        assertEquals(600, tracker.get(yield));
        assertEquals(100.0, tracker.getPercentage(q30, yield), 0.00001);
        TrackerKey undeterminedYield = new TrackerKey(Sample, "UNDETERMINED-Yield");
        assertEquals(50.0, tracker.getPercentage(undeterminedYield, yield), 0.00001);
        TrackerKey laneYield = new TrackerKey(Lane, "L001-Yield");
        assertEquals(600, tracker.get(laneYield));
    }

    @Test
    public void computesCorrectStatsOfUndeterminedLane() throws IOException{
        File dir = new File("./src/test/resources/100UndeterminedBaseCalls1Lane");
        FastqTracker tracker = FastqStats.processDir(dir);
        TrackerKey yield = new TrackerKey(Flowcell, "Yield");
        TrackerKey q30 = new TrackerKey(Flowcell, "Q30");
        assertEquals(600, tracker.get(yield));
        assertEquals(100.0, tracker.getPercentage(q30, yield), 0.00001);
        TrackerKey undeterminedYield = new TrackerKey(Sample, "UNDETERMINED-Yield");
        assertEquals(50.0, tracker.getPercentage(undeterminedYield, yield), 0.00001);
        TrackerKey lane1Yield = new TrackerKey(Lane, "L001-Yield");
        assertEquals(300, tracker.get(lane1Yield));
        TrackerKey lane2Yield = new TrackerKey(Lane, "L002-Yield");
        assertEquals(300, tracker.get(lane2Yield));
        System.out.println(tracker);
    }


}
