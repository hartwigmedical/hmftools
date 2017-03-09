package com.hartwig.hmftools.fastqstats;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class FastqStats {
    private static final Logger LOGGER = LogManager.getLogger(FastqStatsRunner.class);

    public static FastqTracker processFile(String fileName) throws IOException{
        FastqTracker tracker = new FastqTracker();
        File f = new File(fileName);
        FastqData data = processFile(f, tracker);
        String lane = f.getName().split("_")[3];
        tracker.addToFlowcell(data);
        tracker.addToLane(lane, data);
        tracker.addToSample(f.getName(), data);
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
        FastqTracker tracker = new FastqTracker();

        for(File f : files){
            if(f.isDirectory() && f.getName().startsWith("HMFreg")){
                LOGGER.info("Found HMFreg folder: " + f.getName());
                File[] sampleFolders = f.listFiles();
                for(File sampleFolder : sampleFolders){
                    if(sampleFolder.isDirectory()){
                        LOGGER.info("Found sample folder: " + sampleFolder.getName());
                        String sampleName = sampleFolder.getName();
                        for(File fastq: sampleFolder.listFiles()){
                            try {
                                FastqData data = processFile(fastq, tracker);
                                String lane = fastq.getName().split("_")[3];
                                tracker.addToFlowcell(data);
                                tracker.addToLane(lane, data);
                                tracker.addToSample(sampleName, data);
                            }
                            catch (IOException e) {
                                LOGGER.info("Ignoring file: " + fastq.getName());
                            }
                        }
                    }
                }
            }
            else if(!f.isDirectory() && f.getName().startsWith("UNDETERMINED")){
                // undetermined files
                LOGGER.info("Found undetermined file: " + f.getName());
                FastqData data = processFile(f, tracker);
                String lane = f.getName().split("_")[3];
                tracker.addToFlowcell(data);
                tracker.addToUndetermined(data);
                tracker.addToLane(lane, data);
            }
        }
        return tracker;
    }

    private static FastqData processFile(File f, FastqTracker tracker) throws IOException{
        InputStream in;
        int size = 1048576;
        if(f.getName().endsWith(".fastq.gz")) {
            in = new GZIPInputStream(new FileInputStream(new File(f.getCanonicalPath())), size);
        }
        else if(f.getName().endsWith(".fastq")) {
            in = new FileInputStream(new File(f.getCanonicalPath()));
        }
        else {
            throw new IOException("Unrecognized file format.");
        }
        LOGGER.info("Processing file: " + f.getName());
        FastqReader fr = new FastqReader(in, size);
        long startTime = System.currentTimeMillis();
        FastqData data = fr.read();
        long endTime = System.currentTimeMillis();
        fr.close();
        LOGGER.info("Finished processing file: " + f.getName() + " in " + (endTime - startTime) + "ms.");
        return data;
    }
}
