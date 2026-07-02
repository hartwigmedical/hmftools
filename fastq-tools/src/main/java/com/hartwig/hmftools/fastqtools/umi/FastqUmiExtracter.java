package com.hartwig.hmftools.fastqtools.umi;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createGzipBufferedWriter;
import static com.hartwig.hmftools.fastqtools.FastqCommon.APP_NAME;
import static com.hartwig.hmftools.fastqtools.FastqCommon.FQ_LOGGER;
import static com.hartwig.hmftools.fastqtools.FastqCommon.FASTQ_SUFFIX_SHORT;
import static com.hartwig.hmftools.fastqtools.FastqCommon.FASTQ_SUFFIX_STANDARD;
import static com.hartwig.hmftools.fastqtools.FastqCommon.READ_ID_BREAK;
import static com.hartwig.hmftools.fastqtools.FastqCommon.READ_ID_DELIM;
import static com.hartwig.hmftools.fastqtools.FastqCommon.READ_ID_START;
import static com.hartwig.hmftools.fastqtools.FastqCommon.READ_ITEM_BASES;
import static com.hartwig.hmftools.fastqtools.FastqCommon.READ_ITEM_ID;
import static com.hartwig.hmftools.fastqtools.FastqCommon.READ_ITEM_QUALS;
import static com.hartwig.hmftools.fastqtools.FastqCommon.READ_LINE_COUNT;
import static com.hartwig.hmftools.fastqtools.umi.UmiConfig.FASTQ_FILES_DELIM;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class FastqUmiExtracter
{
    private final UmiConfig mConfig;

    private BufferedWriter mWriterR1;
    private BufferedWriter mWriterR2;

    private static final int LINE_LOG_COUNT = 1000000;

    private long mReadCount;
    private final UmiExtractor mUmiExtractor;
    private final KnownUmis mKnownUmis;

    public FastqUmiExtracter(final ConfigBuilder configBuilder)
    {
        mConfig = new UmiConfig(configBuilder);

        mUmiExtractor = new UmiExtractor(mConfig);
        mKnownUmis = new KnownUmis(mConfig);

        mReadCount = 0;
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

        processFiles(fastqFiles[0], fastqFiles[1]);

        mKnownUmis.logResults(mReadCount);

        FQ_LOGGER.info("extraction complete, totalReads({}), mins({})", mReadCount, runTimeMinsStr(startTimeMs));
    }

    private BufferedWriter createOutputWriter(final String inputFile)
    {
        String fastqFile = inputFile.substring(inputFile.lastIndexOf(File.separator) + 1);

        int extensionIndex = fastqFile.contains(FASTQ_SUFFIX_STANDARD) ?
                fastqFile.lastIndexOf("." + FASTQ_SUFFIX_STANDARD) : fastqFile.lastIndexOf("." + FASTQ_SUFFIX_SHORT);

        String outputFile = mConfig.OutputDir + fastqFile.substring(0, extensionIndex) + '.' + mConfig.OutputId + fastqFile.substring(extensionIndex);

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

    private void processFiles(final String r1File, final String r2File)
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
                    if(!processReadBases(r1ReadBuffer, r2ReadBuffer))
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

                    ++mReadCount;
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

    private boolean processReadBases(final String[] r1ReadBuffer, final String[] r2ReadBuffer)
    {
        /* expected format:
        @A00121:853:H5JNNMMABC:1:1101:1443:1047 1:N:0:GGCACAACCT+CAGGAGTCTA
        GNGAGATGGAGAATTTTCTGGAGATGTCTGAGGAATTTTTTCCTCAGTCTTAAGAGTA etc
        +
        F#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFF etc
        */

        // data validation
        if(r1ReadBuffer[READ_ITEM_ID].charAt(0) != READ_ID_START || r2ReadBuffer[READ_ITEM_ID].charAt(0) != READ_ID_START)
            return false;
        
        int minReadLength = max(mUmiExtractor.adapterUmiLength(), mConfig.UmiLength);

        if(r1ReadBuffer[READ_ITEM_BASES].length() <= minReadLength || r1ReadBuffer[READ_ITEM_QUALS].length() <= minReadLength)
            return false;

        int delimIndex = r1ReadBuffer[READ_ITEM_ID].indexOf(READ_ID_BREAK);

        if(delimIndex < 1 || r2ReadBuffer[READ_ITEM_ID].charAt(delimIndex) != READ_ID_BREAK)
            return false;

        String readId1 = r1ReadBuffer[READ_ITEM_ID].substring(1, delimIndex);
        String readId2 = r2ReadBuffer[READ_ITEM_ID].substring(1, delimIndex);

        if(!readId1.equals(readId2))
            return false;

        if(mConfig.AdapterLength > 0)
        {
            mUmiExtractor.adjustWithAdapter(readId1, readId2, delimIndex, r1ReadBuffer,  r2ReadBuffer);
        }
        else if(!mKnownUmis.enabled())
        {
            mKnownUmis.adjustWithKnownUmi(readId1, delimIndex, r1ReadBuffer, r2ReadBuffer);
        }
        else
        {
            mUmiExtractor.adjustWithFixedUmi(readId1, delimIndex, r1ReadBuffer, r2ReadBuffer);
        }

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
