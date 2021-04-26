package com.hartwig.hmftools.fastqstats;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Collection;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.function.Consumer;
import java.util.zip.GZIPInputStream;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.util.concurrent.FutureCallback;
import com.google.common.util.concurrent.Futures;
import com.google.common.util.concurrent.ListenableFuture;
import com.google.common.util.concurrent.ListeningExecutorService;
import com.google.common.util.concurrent.MoreExecutors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class FastqStats {

    private static final Logger LOGGER = LogManager.getLogger(FastqStats.class);

    /**
     * Counts yield and q30 of fastqs in the fastqsPerSample multimap, using 1 thread per file.
     * The yield and q30 of the Undetermined sample will count towards the total yield and q30 of the flowcell.
     *
     * @param fastqsPerSample multimap of sampleName and fastqs to process
     * @param threadCount     number of maximum threads
     * @return FastqTracker with yield and q30 stats for the fastqs processed.
     */

    @NotNull
    static FastqTracker processFastqs(@NotNull final Multimap<String, File> fastqsPerSample, final int threadCount)
            throws InterruptedException {
        LOGGER.info("Using {} threads. Processing {} fastQ files.", threadCount, fastqsPerSample.size());
        final FastqTrackerWrapper tracker = new FastqTrackerWrapper();
        final ListeningExecutorService threadPool = MoreExecutors.listeningDecorator(Executors.newFixedThreadPool(threadCount));

        for (final String sampleName : fastqsPerSample.keySet()) {
            final Collection<File> fastqs = fastqsPerSample.get(sampleName);
            for (final File fastq : fastqs) {
                final String laneName = getLaneName(fastq);
                final ListenableFuture<FastqData> futureResult = threadPool.submit(() -> processFile(fastq));
                addCallback(futureResult, (data) -> tracker.addDataFromSampleFile(sampleName, laneName, data),
                        (error) -> LOGGER.error("Failed to process file: {}", fastq.getName(), error));
            }
        }
        threadPool.shutdown();
        threadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        return tracker.tracker();
    }

    /**
     * Looks for folders starting with HMFreg in the current baseCallsDir. If any folder is found, assumes all subfolders represent
     * samples and collects all Fastq files found in the samples subfolders.
     * Any Fastq files in the baseCallsDir that start with 'Undetermined' will be assigned to an 'Undetermined' sample.
     *
     * @param baseCallsDir BaseCalls directory
     * @return a multimap with sample names as keys and fastq files as values.
     * @throws IOException if an error occurs when reading HMFReg/sample folders or fastq files.
     */

    @NotNull
    static Multimap<String, File> getFastqsFromBaseCallsDir(@NotNull final File baseCallsDir) throws IOException {
        final Multimap<String, File> fastqsPerSample = ArrayListMultimap.create();
        final File[] files = baseCallsDir.listFiles();
        if (files == null) {
            throw new IOException("List files in " + baseCallsDir.getName() + " returned null.");
        }
        for (final File file : files) {
            if (file.isDirectory() && file.getName().startsWith("HMFreg")) {
                LOGGER.info("Found HMFreg folder: {}", file.getName());
                final File[] sampleFolders = file.listFiles(File::isDirectory);
                if (sampleFolders == null) {
                    throw new IOException("List folders in " + file.getName() + " returned null.");
                }
                for (final File sampleFolder : sampleFolders) {
                    LOGGER.info("Found sample folder: {}", sampleFolder.getName());
                    final File[] fastqFiles =
                            sampleFolder.listFiles((parentDir, fileName) -> fileName.endsWith(".fastq.gz") || fileName.endsWith(".fastq"));
                    if (fastqFiles == null) {
                        throw new IOException("List fastq files in " + sampleFolder.getName() + " returned null.");
                    }
                    fastqsPerSample.putAll(sampleFolder.getName(), Lists.newArrayList(fastqFiles));
                }
            } else if (!file.isDirectory() && file.getName().startsWith("Undetermined") && (file.getName().endsWith(".fastq.gz")
                    || file.getName().endsWith(".fastq"))) {
                LOGGER.info("Found undetermined file: {}", file.getName());
                fastqsPerSample.put("Undetermined", file);
            }
        }
        return fastqsPerSample;
    }

    @NotNull
    static Multimap<String, File> getSingleFastq(@NotNull final String filePath) {
        final Multimap<String, File> fastqsPerSample = ArrayListMultimap.create();
        final File file = new File(filePath);
        fastqsPerSample.put(file.getName(), file);
        return fastqsPerSample;
    }

    @NotNull
    static Multimap<String, File> getFastqsFromDir(@NotNull final File fastqDir) throws IOException {
        final Multimap<String, File> fastqsPerSample = ArrayListMultimap.create();
        final File[] fastqFiles = fastqDir.listFiles(file -> file.getName().endsWith(".fastq.gz") || file.getName().endsWith(".fastq"));
        if (fastqFiles == null) {
            throw new IOException("List files in " + fastqDir.getName() + " returned null.");
        }
        for (final File fastq : fastqFiles) {
            fastqsPerSample.put(fastq.getName(), fastq);
        }
        return fastqsPerSample;
    }

    @NotNull
    private static FastqData processFile(@NotNull final File file) throws IOException {
        final int size = 1048576;
        final InputStream inputStream;
        if (file.getName().endsWith(".fastq.gz")) {
            inputStream = new GZIPInputStream(new FileInputStream(new File(file.getCanonicalPath())), size);
        } else if (file.getName().endsWith(".fastq")) {
            inputStream = new FileInputStream(new File(file.getCanonicalPath()));
        } else {
            throw new IOException("Unrecognized file format.");
        }
        LOGGER.info("Processing file: {}", file.getName());
        final FastqReader fastqReader = new FastqReader(inputStream, size);
        final long startTime = System.currentTimeMillis();
        final FastqData data = fastqReader.read();
        final long endTime = System.currentTimeMillis();
        fastqReader.close();

        LOGGER.info("Finished processing file: {} in {}ms.", file.getName(), endTime - startTime);
        return data;
    }

    private static <T> void addCallback(@NotNull final ListenableFuture<T> future, @NotNull final Consumer<T> onSuccess,
            @NotNull final Consumer<Throwable> onFailure) {
        Futures.addCallback(future, new FutureCallback<T>() {
            @Override
            public void onSuccess(@Nullable final T t) {
                assert t != null;
                onSuccess.accept(t);
            }

            @Override
            public void onFailure(@NotNull final Throwable throwable) {
                onFailure.accept(throwable);
            }
        }, MoreExecutors.directExecutor());
    }

    /**
     * Get the name of the lane from the file name.
     * Assumes the name of the file has the following, underscore-separated format: STUDYCODE_FLOWCELL_S1_L001_R1_001.fastq
     * where the lane name is on the 4th position.
     *
     * @param file file to extract the lane name from
     * @return lane name, or "Unknown" if the lane name could not be extracted
     */
    @NotNull
    private static String getLaneName(@NotNull File file) {
        final String[] fileNameArray = file.getName().split("_");
        if (fileNameArray.length >= 4) {
            return fileNameArray[3];
        } else {
            LOGGER.warn("Could not get lane name from {}", file.getName());
            return "Unknown";
        }
    }
}
