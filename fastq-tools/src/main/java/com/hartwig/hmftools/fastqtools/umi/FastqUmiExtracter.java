package com.hartwig.hmftools.fastqtools.umi;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createGzipBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
import static com.hartwig.hmftools.fastqtools.FastqCommon.APP_NAME;
import static com.hartwig.hmftools.fastqtools.FastqCommon.FASTQ_ZIP_EXTENSION;
import static com.hartwig.hmftools.fastqtools.FastqCommon.FQ_LOGGER;
import static com.hartwig.hmftools.fastqtools.FastqCommon.READ_ITEM_QUALS;
import static com.hartwig.hmftools.fastqtools.FastqCommon.READ_LINE_COUNT;
import static com.hartwig.hmftools.fastqtools.umi.UmiConfig.FASTQ_FILES_DELIM;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class FastqUmiExtracter
{
    private final UmiConfig mConfig;

    private BufferedWriter mWriterR1;
    private BufferedWriter mWriterR2;

    private static final int LINE_LOG_COUNT = 1000000;

    private long mLineCount;
    private final UmiExtractor mUmiExtractor;

    public FastqUmiExtracter(final ConfigBuilder configBuilder)
    {
        mConfig = new UmiConfig(configBuilder);

        mUmiExtractor = new UmiExtractor(mConfig);

        mLineCount = 0;
    }

    public void run()
    {
        if(mConfig.OutputDir == null || mConfig.FastqFiles == null)
            System.exit(1);

        String[] fastqFiles = null;
        if(mConfig.FastqFiles.contains(FASTQ_FILES_DELIM))
        {
            fastqFiles = mConfig.FastqFiles.split(FASTQ_FILES_DELIM, 2);
        }
        else if(mConfig.FastqFiles.contains("*"))
        {
            fastqFiles = new String[2];
            fastqFiles[0] = mConfig.FastqFiles.replaceAll("\\*", "1");
            fastqFiles[1] = mConfig.FastqFiles.replaceAll("\\*", "2");
        }

        if(fastqFiles.length != 2)
        {
            FQ_LOGGER.info("invalid fastq file config: {}", mConfig.FastqFiles);
            System.exit(1);
        }

        if(!Files.exists(Paths.get(fastqFiles[0])) || !Files.exists(Paths.get(fastqFiles[1])))
        {
            FQ_LOGGER.info("fastq files do not exist: {}", mConfig.FastqFiles);
            System.exit(1);
        }

        FQ_LOGGER.info("starting Fastq UMI extractions with files: {}", mConfig.FastqFiles);

        long startTimeMs = System.currentTimeMillis();

        try
        {
            processFiles(fastqFiles[0], fastqFiles[1]);
        }
        catch(Exception e)
        {
            FQ_LOGGER.error("file processing failed: {}", e.toString());
            e.printStackTrace();
        }

        long readCount = mLineCount / 4;
        mUmiExtractor.logResults(readCount);

        FQ_LOGGER.info("extraction complete, totalReads({}), mins({})", readCount, runTimeMinsStr(startTimeMs));
    }

    private BufferedWriter createOutputWriter(final String inputFile)
    {
        String outputFile = mConfig.outputFastqFilename(inputFile);

        try
        {
            return outputFile.endsWith(".gz") ? createGzipBufferedWriter(outputFile) : createBufferedWriter(outputFile);
        }
        catch(IOException e)
        {
            FQ_LOGGER.error("error creating fastq output file({}): {}", outputFile, e.toString());
            System.exit(1);
            return null;
        }
    }

    private static BufferedReader openReader(final String fastqFile) throws Exception
    {
        InputStream fis = Files.newInputStream(Paths.get(fastqFile));
        BufferedInputStream bis = new BufferedInputStream(fis);
        GZIPInputStream gzis = new GZIPInputStream(bis);
        InputStreamReader isr = new InputStreamReader(gzis, StandardCharsets.UTF_8);
        return new BufferedReader(isr);
    }

    private void processFiles(final String inputFileR1, final String inputFileR2) throws Exception
    {
        ExecutorService executor = Executors.newFixedThreadPool(mConfig.Threads);
        List<Future<Void>> futures = new ArrayList<>();

        String outputFileR1 = mConfig.outputFastqFilename(inputFileR1);
        String outputFileR2 = mConfig.outputFastqFilename(inputFileR2);
        String outputPrefixR1 = outputFileR1.substring(0, outputFileR1.lastIndexOf(FASTQ_ZIP_EXTENSION));
        String outputPrefixR2 = outputFileR2.substring(0, outputFileR2.lastIndexOf(FASTQ_ZIP_EXTENSION));

        int linesPerChunk = mConfig.ReadsPerChunk * READ_LINE_COUNT;

        BufferedReader readerR1 = openReader(inputFileR1);
        BufferedReader readerR2 = openReader(inputFileR2);

        String[] r1Lines = new String[linesPerChunk];
        String[] r2Lines = new String[linesPerChunk];

        String line;
        int chunkIndex = 0;
        int lineIndex = 0;

        int nextLogCount = LINE_LOG_COUNT;

        while((line = readerR1.readLine()) != null)
        {
            r1Lines[lineIndex] = line;
            r2Lines[lineIndex] = readerR2.readLine(); // assumes equal lines
            ++lineIndex;
            ++mLineCount;

            if(mLineCount >= nextLogCount)
            {
                nextLogCount += LINE_LOG_COUNT;
                FQ_LOGGER.debug("processed {} lines", mLineCount);
            }

            // accumulate all lines for a single read
            if(lineIndex == linesPerChunk)
            {
                int currentIndex = chunkIndex++;

                Path outputTempFileR1 = formTemporaryFilename(outputPrefixR1, currentIndex);
                Path outputTempFileR2 = formTemporaryFilename(outputPrefixR2, currentIndex);

                String[] r1LinesCp = r1Lines;
                String[] r2LinesCp = r2Lines;

                Future<Void> task = executor.submit(() ->
                {
                    processReadGroupLines(r1LinesCp, r2LinesCp);
                    compressTextLines(r1LinesCp, outputTempFileR1);
                    compressTextLines(r2LinesCp, outputTempFileR2);
                    return null;
                });

                futures.add(task);

                // reset batch array for the next file chunk
                lineIndex = 0;

                // clear previous data since buffers are passed to threads
                r1Lines = new String[linesPerChunk];
                r2Lines = new String[linesPerChunk];
            }
        }

        // handle remaining reads (less than chunk size
        if(lineIndex > 0)
        {
            Path outputTempFileR1 = formTemporaryFilename(outputPrefixR1, chunkIndex);
            Path outputTempFileR2 = formTemporaryFilename(outputPrefixR2, chunkIndex);

            String[] r1LinesCp = r1Lines;
            String[] r2LinesCp = r2Lines;

            Future<Void> task = executor.submit(() ->
            {
                processReadGroupLines(r1LinesCp, r2LinesCp);
                compressTextLines(r1LinesCp, outputTempFileR1);
                compressTextLines(r2LinesCp, outputTempFileR2);
                return null;
            });

            futures.add(task);
        }

        for(Future<Void> future : futures)
        {
            future.get();
        }

        executor.shutdown();

        FQ_LOGGER.info("read group processing complete");

        if(mConfig.Threads > 1)
        {
            FQ_LOGGER.info("merging {} temporary files", mConfig.Threads);
        }

        mergeGzipBlocks(outputPrefixR1, outputFileR1);
        mergeGzipBlocks(outputPrefixR2, outputFileR2);
    }

    private static Path formTemporaryFilename(final String outputPrefix, int fileIndex)
    {
        return Paths.get(format("%s_%05d.gz", outputPrefix, fileIndex));
    }

    private void processReadGroupLines(final String[] r1Lines, final String[] r2Lines)
    {
        String[] r1ReadBuffer = new String[READ_ITEM_QUALS+1];
        String[] r2ReadBuffer = new String[READ_ITEM_QUALS+1];

        for(int i = 0; i < r1Lines.length; i += READ_LINE_COUNT)
        {
            for(int j = 0; j < READ_LINE_COUNT; ++j)
            {
                if(r1Lines[i + j] == null)
                    return;

                r1ReadBuffer[j] = r1Lines[i + j];
                r2ReadBuffer[j] = r2Lines[i + j];
            }

            mUmiExtractor.processReadBases(r1ReadBuffer, r2ReadBuffer);

            // copy back
            for(int j = 0; j < READ_LINE_COUNT; ++j)
            {
                r1Lines[i + j] = r1ReadBuffer[j];
                r2Lines[i + j] = r2ReadBuffer[j];
            }
        }
    }

    private static void compressTextLines(final String[] lines, final Path outputPath) throws IOException
    {
        try(OutputStream fos = Files.newOutputStream(outputPath);
            BufferedOutputStream bos = new BufferedOutputStream(fos);
            GZIPOutputStream gzos = new GZIPOutputStream(bos);
            BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(gzos, StandardCharsets.UTF_8)))
        {
            for(String line : lines)
            {
                if(line == null)
                    break;

                writer.write(line);
                writer.newLine();
            }
        }
    }

    private static final String ZIPPED_EXTENSION = ".gz";

    private void mergeGzipBlocks(final String pathPrefix, final String combinedOutputFile) throws IOException
    {
        String filePrefix = filenamePart(pathPrefix) + "_";

        // scan the directory for target parts and sort them numerically by sequence ID
        List<Path> interimFiles = Files.list(Paths.get(mConfig.OutputDir))
                .filter(p -> p.getFileName().toString().startsWith(filePrefix) && p.toString().endsWith(ZIPPED_EXTENSION))
                .collect(Collectors.toList());

        Collections.sort(interimFiles);

        if(interimFiles.isEmpty())
        {
            throw new FileNotFoundException("No split chunks found matching the pattern: " + pathPrefix + "_*.gz");
        }

        // merge raw compressed byte segments together sequentially
        OutputStream destStream = Files.newOutputStream(Paths.get(combinedOutputFile));
        BufferedOutputStream bos = new BufferedOutputStream(destStream);

        byte[] buffer = new byte[8192];

        for(Path part : interimFiles)
        {
            InputStream sourceStream = Files.newInputStream(part);
            BufferedInputStream bis = new BufferedInputStream(sourceStream);

            int bytesRead;

            while((bytesRead = bis.read(buffer)) != -1)
            {
                bos.write(buffer, 0, bytesRead);
            }
        }

        bos.flush();

        // clean-up interim files
        for(Path interimFile : interimFiles)
        {
            Files.delete(interimFile);
        }
    }
    private void processFilesOld(final String r1File, final String r2File)
    {
        try
        {
            mWriterR1 = createOutputWriter(r1File);
            mWriterR2 = createOutputWriter(r2File);

            BufferedReader r1Reader = createBufferedReader(r1File);
            BufferedReader r2Reader = createBufferedReader(r2File);

            int lineCount = 0;

            int readLineCount = 0;
            String[] r1ReadBuffer = new String[READ_ITEM_QUALS+1];
            String[] r2ReadBuffer = new String[READ_ITEM_QUALS+1];

            while(true)
            {
                r1ReadBuffer[readLineCount] = r1Reader.readLine();
                r2ReadBuffer[readLineCount] = r2Reader.readLine();

                if(r1ReadBuffer[readLineCount] == null || r2ReadBuffer[readLineCount] == null)
                    break;

                ++readLineCount;

                if(readLineCount == READ_LINE_COUNT)
                {
                    if(!processReadBasesOld(r1ReadBuffer, r2ReadBuffer))
                    {
                        FQ_LOGGER.error("invalid entries at line({})", lineCount);

                        for(int i = 0; i < r1ReadBuffer.length; ++i)
                        {
                            FQ_LOGGER.error("R1 item {}: {}", i, r1ReadBuffer[i]);
                            FQ_LOGGER.error("R2 item {}: {}", i, r2ReadBuffer[i]);
                        }

                        System.exit(1);
                    }

                    readLineCount = 0;

                    ++mLineCount;
                }

                ++lineCount;

                if(lineCount > 0 && (lineCount % LINE_LOG_COUNT) == 0)
                {
                    FQ_LOGGER.info("processed {} lines", lineCount);
                }
            }

            mWriterR1.close();
            mWriterR2.close();
        }
        catch(IOException e)
        {
            FQ_LOGGER.error("error reading fastq({}): {}", mConfig.FastqFiles, e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }

    private boolean processReadBasesOld(final String[] r1ReadBuffer, final String[] r2ReadBuffer)
    {
        mUmiExtractor.processReadBases(r1ReadBuffer, r2ReadBuffer);

        try
        {
            for(int i = 0; i < r1ReadBuffer.length; ++i)
            {
                mWriterR1.write(r1ReadBuffer[i]);
                mWriterR1.newLine();

                mWriterR2.write(r2ReadBuffer[i]);
                mWriterR2.newLine();
            }

            return true;
        }
        catch(IOException e)
        {
            FQ_LOGGER.error("failed to write output file: {}", e.toString());
            return false;
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        UmiConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        FastqUmiExtracter fastqUmiExtracter = new FastqUmiExtracter(configBuilder);
        fastqUmiExtracter.run();
    }
}
