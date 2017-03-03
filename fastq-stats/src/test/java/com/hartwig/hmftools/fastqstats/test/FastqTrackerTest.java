package com.hartwig.hmftools.fastqstats.test;

import static com.hartwig.hmftools.fastqstats.TrackerType.Flowcell;
import org.junit.Test;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.TreeMap;
import static org.junit.Assert.assertEquals;
import com.hartwig.hmftools.fastqstats.FastqReader;
import com.hartwig.hmftools.fastqstats.PredicateTracker;
import com.hartwig.hmftools.fastqstats.FastqTracker;
import com.hartwig.hmftools.fastqstats.Tracker;
import com.hartwig.hmftools.fastqstats.TrackerKey;
import com.hartwig.hmftools.fastqstats.YieldTracker;

public class FastqTrackerTest {
    @Test
    public void computesCorrectPercentageFor100PercentFile() throws IOException{
        FileInputStream fin = new FileInputStream(new File("./src/test/resources/fastq/q30-100.fastq"));
        FastqTracker tracker = createTestTracker();
        FastqReader fr = new FastqReader(fin, tracker);
        TrackerKey yield = new TrackerKey(Flowcell, "Yield");
        TrackerKey q30 = new TrackerKey(Flowcell, "Q30");
        fr.read(new TrackerKey[]{yield, q30});
        double q30Percentage = tracker.getPercentage(q30, yield);
        assertEquals(100, tracker.get(yield));
        assertEquals(100.0, q30Percentage, 0.00001);
    }

    @Test
    public void computesCorrectPercentageFor0PercentFile() throws IOException{
        FileInputStream fin = new FileInputStream(new File("./src/test/resources/fastq/q30-0.fastq"));
        FastqTracker tracker = createTestTracker();
        FastqReader fr = new FastqReader(fin, tracker);
        TrackerKey yield = new TrackerKey(Flowcell, "Yield");
        TrackerKey q30 = new TrackerKey(Flowcell, "Q30");
        fr.read(new TrackerKey[]{yield, q30});
        double q30Percentage = tracker.getPercentage(q30, yield);
        assertEquals(100, tracker.get(yield));
        assertEquals(0.0, q30Percentage,0.00001);
    }

    @Test
    public void computesCorrectPercentageFor10PercentFile() throws IOException{
        FileInputStream fin = new FileInputStream(new File("./src/test/resources/fastq/q30-10.fastq"));
        FastqTracker tracker = createTestTracker();
        FastqReader fr = new FastqReader(fin, tracker);
        TrackerKey yield = new TrackerKey(Flowcell, "Yield");
        TrackerKey q30 = new TrackerKey(Flowcell, "Q30");
        fr.read(new TrackerKey[]{yield, q30});
        double q30Percentage = tracker.getPercentage(q30, yield);
        assertEquals(100, tracker.get(yield));
        assertEquals(10.0, q30Percentage, 0.00001);
    }

    private FastqTracker createTestTracker(){
        return new FastqTracker(new TreeMap<TrackerKey, Tracker>(){{
            put(new TrackerKey(Flowcell, "Yield"), new YieldTracker());
            put(new TrackerKey(Flowcell, "Q30"), new PredicateTracker(x -> x >= 30));
        }});
    }
}
