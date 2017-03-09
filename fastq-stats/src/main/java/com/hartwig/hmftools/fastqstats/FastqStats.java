package com.hartwig.hmftools.fastqstats;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class FastqStats {
    private static final Logger LOGGER = LogManager.getLogger(FastqStatsRunner.class);

    @NotNull
    public static FastqTracker processFile(@NotNull String filePath) throws IOException {
        FastqTracker tracker = new FastqTracker();
        File file = new File(filePath);
        FastqData data = processFile(file);
        String lane = file.getName().split("_")[3];
        tracker.addToFlowcell(data);
        tracker.addToLane(lane, data);
        tracker.addToSample(file.getName(), data);
        return tracker;
    }

    /**
     * Looks for folders starting with HMFreg in the current (BaseCalls) dir. If any folder with
     * that pattern is found, assumes all subfolders represent samples and processes all Fastq
     * files found in these subfolders.
     * Any Fastq files found in the current (BaseCalls) dir that have the string UNDETERMINED at
     * the start will be assigned to an 'Undetermined' sample. The yield and q30 from these
     * files will count towards the total yield and q30 of the flowcell.
     *
     * @param dir BaseCalls directory
     */
    @NotNull
    public static FastqTracker processDir(@NotNull File dir) throws IOException {
        File[] files = dir.listFiles();
        if (files == null) {
            throw new IOException("List files in " + dir.getName() + " returned null.");
        }
        FastqTracker tracker = new FastqTracker();

        for (final File file : files) {
            if (file.isDirectory() && file.getName().startsWith("HMFreg")) {
                LOGGER.info("Found HMFreg folder: " + file.getName());
                File[] sampleFolders = file.listFiles(File::isDirectory);
                if (sampleFolders == null) {
                    throw new IOException("List folders in " + file.getName() + " returned null.");
                }
                for (final File sampleFolder : sampleFolders) {
                    LOGGER.info("Found sample folder: " + sampleFolder.getName());
                    String sampleName = sampleFolder.getName();
                    File[] fastqFiles = sampleFolder.listFiles(
                            (parentDir, fileName) -> fileName.endsWith(".fastq.gz") || fileName.endsWith(".fastq"));
                    if (fastqFiles == null) {
                        throw new IOException("List fastq files in " + sampleFolder.getName() + " returned null.");
                    }
                    for (final File fastq : fastqFiles) {
                        FastqData data = processFile(fastq);
                        String lane = fastq.getName().split("_")[3];
                        tracker.addToFlowcell(data);
                        tracker.addToLane(lane, data);
                        tracker.addToSample(sampleName, data);
                    }
                }
            } else if (!file.isDirectory() && file.getName().startsWith("UNDETERMINED") && (
                    file.getName().endsWith(".fastq.gz") || file.getName().endsWith(".fastq"))) {
                // undetermined files
                LOGGER.info("Found undetermined file: " + file.getName());
                FastqData data = processFile(file);
                String lane = file.getName().split("_")[3];
                tracker.addToFlowcell(data);
                tracker.addToUndetermined(data);
                tracker.addToLane(lane, data);
            }
        }
        return tracker;
    }

    @NotNull
    private static FastqData processFile(@NotNull File f) throws IOException {
        int size = 1048576;
        final InputStream in;
        if (f.getName().endsWith(".fastq.gz")) {
            in = new GZIPInputStream(new FileInputStream(new File(f.getCanonicalPath())), size);
        } else if (f.getName().endsWith(".fastq")) {
            in = new FileInputStream(new File(f.getCanonicalPath()));
        } else {
            throw new IOException("Unrecognized file format.");
        }
        LOGGER.info("Processing file: " + f.getName());
        final FastqReader fr = new FastqReader(in, size);
        long startTime = System.currentTimeMillis();
        final FastqData data = fr.read();
        long endTime = System.currentTimeMillis();
        fr.close();
        LOGGER.info("Finished processing file: " + f.getName() + " in " + (endTime - startTime) + "ms.");
        return data;
    }
}
