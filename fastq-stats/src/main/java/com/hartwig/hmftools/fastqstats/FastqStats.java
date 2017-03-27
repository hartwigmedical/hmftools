package com.hartwig.hmftools.fastqstats;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.function.Consumer;
import java.util.zip.GZIPInputStream;

import com.google.common.util.concurrent.FutureCallback;
import com.google.common.util.concurrent.Futures;
import com.google.common.util.concurrent.ListenableFuture;
import com.google.common.util.concurrent.ListeningExecutorService;
import com.google.common.util.concurrent.MoreExecutors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class FastqStats {
    private static final Logger LOGGER = LogManager.getLogger(FastqStatsRunner.class);

    @NotNull
    static FastqTracker processFile(@NotNull String filePath) throws IOException {
        final FastqTracker tracker = new FastqTracker();
        final File file = new File(filePath);
        final FastqData data = processFile(file);
        final String lane = file.getName().split("_")[3];
        return tracker.addToFlowcell(data).addToSample(file.getName(), lane, data);
    }

    /**
     * Looks for folders starting with HMFreg in the current (BaseCalls) dir. If any folder with
     * that pattern is found, assumes all subfolders represent samples and processes all Fastq
     * files found in these subfolders.
     * Any Fastq files found in the current (BaseCalls) dir that have the string UNDETERMINED at
     * the start will be assigned to an 'Undetermined' sample. The yield and q30 from these
     * files will count towards the total yield and q30 of the flowcell.
     *
     * @param dir         BaseCalls directory
     * @param threadCount number of maximum threads to use
     * @return
     * @throws IOException
     * @throws InterruptedException
     */
    @NotNull
    static FastqTracker processDir(@NotNull File dir, int threadCount) throws IOException, InterruptedException {
        final File[] files = dir.listFiles();
        if (files == null) {
            throw new IOException("List files in " + dir.getName() + " returned null.");
        }
        final FastqTrackerWrapper tracker = new FastqTrackerWrapper();
        LOGGER.info("Using " + threadCount + " threads.");
        final ListeningExecutorService threadPool = MoreExecutors.listeningDecorator(
                Executors.newFixedThreadPool(threadCount));

        for (final File file : files) {
            if (file.isDirectory() && file.getName().startsWith("HMFreg")) {
                LOGGER.info("Found HMFreg folder: " + file.getName());
                final File[] sampleFolders = file.listFiles(File::isDirectory);
                if (sampleFolders == null) {
                    throw new IOException("List folders in " + file.getName() + " returned null.");
                }
                for (final File sampleFolder : sampleFolders) {
                    LOGGER.info("Found sample folder: " + sampleFolder.getName());
                    final File[] fastqFiles = sampleFolder.listFiles(
                            (parentDir, fileName) -> fileName.endsWith(".fastq.gz") || fileName.endsWith(".fastq"));
                    if (fastqFiles == null) {
                        throw new IOException("List fastq files in " + sampleFolder.getName() + " returned null.");
                    }
                    for (final File fastq : fastqFiles) {
                        final String lane = fastq.getName().split("_")[3];
                        final ListenableFuture<FastqData> futureResult = threadPool.submit(() -> processFile(fastq));
                        addCallback(futureResult,
                                (data) -> tracker.addDataFromSampleFile(sampleFolder.getName(), lane, data),
                                (error) -> LOGGER.error("Failed to process file: " + fastq.getName(), error));
                    }
                }
            } else if (!file.isDirectory() && file.getName().startsWith("Undetermined") && (
                    file.getName().endsWith(".fastq.gz") || file.getName().endsWith(".fastq"))) {
                LOGGER.info("Found undetermined file: " + file.getName());
                final String lane = file.getName().split("_")[3];
                final ListenableFuture<FastqData> futureResult = threadPool.submit(() -> processFile(file));
                addCallback(futureResult, (data) -> tracker.addDataFromSampleFile("Undetermined", lane, data),
                        (error) -> LOGGER.error("Failed to process file: " + file.getName(), error));
            }
        }
        threadPool.shutdown();
        threadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        return tracker.tracker();
    }

    @NotNull
    private static FastqData processFile(@NotNull File f) throws IOException {
        final int size = 1048576;
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
        final long startTime = System.currentTimeMillis();
        final FastqData data = fr.read();
        final long endTime = System.currentTimeMillis();
        fr.close();
        LOGGER.info("Finished processing file: " + f.getName() + " in " + (endTime - startTime) + "ms.");
        return data;
    }

    private static <T> void addCallback(@NotNull ListenableFuture<T> future, @NotNull Consumer<T> onSuccess,
            @NotNull Consumer<Throwable> onFailure) {
        Futures.addCallback(future, new FutureCallback<T>() {
            @Override
            public void onSuccess(@NotNull final T t) {
                onSuccess.accept(t);
            }

            @Override
            public void onFailure(@NotNull final Throwable throwable) {
                onFailure.accept(throwable);
            }
        });
    }
}
