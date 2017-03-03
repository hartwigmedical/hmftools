package com.hartwig.hmftools.fastqstats;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;
import static com.hartwig.hmftools.fastqstats.TrackerType.*;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class FastqStats {
    private static final Logger LOGGER = LogManager.getLogger(FastqStatsRunner.class);

    public static FastqTracker processFile(String fileName) throws IOException{
        FastqTracker tracker = createDefaultFastqTracker();
        tracker.putIfAbsent(new TrackerKey(Sample, fileName + "-Yield"), new YieldTracker());
        tracker.putIfAbsent(new TrackerKey(Sample, fileName + "-Q30"), new PredicateTracker(x -> x >= 30));
        processFile(fileName, new File(fileName), tracker);
        return tracker;
    }

    /**
     * Looks for folders starting with HMFreg in the current (BaseCalls) dir. If any folder with
     *  that pattern is found, assumes all subfolders represent samples and processes all Fastq
     *  files found in these subfolders.
     * Any Fastq files found in the current (BaseCalls) dir that have the string UNDETERMINED at
     *  the start will be assigned to an 'Undetermined' sample. The yield and q30 from these
     *  files will count towards the total yield and q30 of the flowcell.
     *
     * @param dir BaseCalls directory
     */
    public static FastqTracker processDir(File dir) throws IOException{
        File[] files = dir.listFiles();
        FastqTracker tracker = createDefaultFastqTracker();

        for(File f : files){
            if(f.isDirectory() && f.getName().startsWith("HMFreg")){
                LOGGER.info("Found HMFreg folder: " + f.getName());
                File[] sampleFolders = f.listFiles();
                for(File sampleFolder : sampleFolders){
                    if(sampleFolder.isDirectory()){
                        LOGGER.info("Found sample folder: " + sampleFolder.getName());
                        String sampleName = sampleFolder.getName();
                        tracker.putIfAbsent(new TrackerKey(Sample, sampleName + "-Yield"), new YieldTracker());
                        tracker.putIfAbsent(new TrackerKey(Sample, sampleName + "-Q30"), new PredicateTracker(x -> x >= 30));
                        for(File fastq: sampleFolder.listFiles()){
                            processFile(sampleName, fastq, tracker);
                        }
                    }
                }
            }
            else if(!f.isDirectory() && f.getName().startsWith("UNDETERMINED")){
                // undetermined files
                LOGGER.info("Found undetermined file: " + f.getName());
                String sampleName = "UNDETERMINED";
                tracker.putIfAbsent(new TrackerKey(Sample, "UNDETERMINED-Yield"), new YieldTracker());
                tracker.putIfAbsent(new TrackerKey(Sample, "UNDETERMINED-Q30"), new PredicateTracker(x -> x >= 30));
                processFile(sampleName, f, tracker);
            }
        }
        return tracker;
    }

    private static void processFile(String sample, File f, FastqTracker tracker) throws IOException{
        InputStream in;
        if(f.getName().endsWith(".fastq.gz")) {
            in = new GZIPInputStream(new FileInputStream(new File(f.getCanonicalPath())));
        }
        else if(f.getName().endsWith(".fastq")) {
            in = new FileInputStream(new File(f.getCanonicalPath()));
        }
        else {
            return;
        }
        LOGGER.info("Processing file: " + f.getName());
        String lane = f.getName().split("_")[3];
        tracker.putIfAbsent(new TrackerKey(Lane, lane + "-Yield"), new YieldTracker());
        tracker.putIfAbsent(new TrackerKey(Lane, lane + "-Q30"), new PredicateTracker(x -> x >= 30));
        FastqReader fr = new FastqReader(in, tracker);
        TrackerKey[] keys = new TrackerKey[] {
            new TrackerKey(Flowcell, "Yield"),
            new TrackerKey(Flowcell, "Q30"),
            new TrackerKey(Lane, lane+"-Yield"),
            new TrackerKey(Lane, lane+"-Q30"),
            new TrackerKey(Sample, sample+"-Yield"),
            new TrackerKey(Sample, sample+"-Q30")
        };
        fr.read(keys);
        fr.close();
        LOGGER.info("Finished processing file: " + f.getName());
    }

    private static FastqTracker createDefaultFastqTracker(){
        return new FastqTracker(new TreeMap<TrackerKey, Tracker>(){{
            put(new TrackerKey(Flowcell, "Yield"), new YieldTracker());
            put(new TrackerKey(Flowcell, "Q30"), new PredicateTracker(x -> x >= 30));
        }});
    }
}
