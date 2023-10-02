package com.hartwig.hmftools.markdups.fastq;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createGzipBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.APP_NAME;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.common.Constants.DEFAULT_DUPLEX_UMI_DELIM;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigItemType;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import org.jetbrains.annotations.NotNull;

public class FastqUmiExtracter
{
    private final String mFastqFiles;
    private final String mOutputDir;
    private final String mOutputId;
    private final int mUmiLength;

    private BufferedWriter mWriterR1;
    private BufferedWriter mWriterR2;

    // config
    private static final String FASTQ_FILES = "fastq_files";
    private static final String UMI_LENGTH = "umi_length";

    private static final String FASTQ_FILES_DELIM = ";";
    private static final int LINE_LOG_COUNT = 1000000;

    public FastqUmiExtracter(final ConfigBuilder configBuilder)
    {
        mFastqFiles = configBuilder.getValue(FASTQ_FILES);
        mOutputDir = parseOutputDir(configBuilder);
        mOutputId = configBuilder.getValue(OUTPUT_ID);
        mUmiLength = configBuilder.getInteger(UMI_LENGTH);
    }

    public void run()
    {
        if(mOutputDir == null || mFastqFiles == null)
            System.exit(1);

        String[] fastqFiles = null;
        if(mFastqFiles.contains(FASTQ_FILES_DELIM))
        {
            fastqFiles = mFastqFiles.split(FASTQ_FILES_DELIM, 2);
        }
        else if(mFastqFiles.contains("*"))
        {
            fastqFiles = new String[2];
            fastqFiles[0] = mFastqFiles.replaceAll("\\*", "1");
            fastqFiles[1] = mFastqFiles.replaceAll("\\*", "2");
        }

        if(fastqFiles.length != 2)
        {
            MD_LOGGER.info("invalid fastq file config: {}", mFastqFiles);
            System.exit(1);
        }

        if(!Files.exists(Paths.get(fastqFiles[0])) || !Files.exists(Paths.get(fastqFiles[1])))
        {
            MD_LOGGER.info("fastq files do not exist: {}", mFastqFiles);
            System.exit(1);
        }

        MD_LOGGER.info("Starting Fastq UMI extractions with files: {}", mFastqFiles);

        long startTimeMs = System.currentTimeMillis();

        processFiles(fastqFiles[0], fastqFiles[1]);

        MD_LOGGER.info("Extraction complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private static final int READ_ITEM_ID = 0;
    private static final int READ_ITEM_BASES = 1;
    private static final int READ_ITEM_SPARE = 2;
    private static final int READ_ITEM_QUALS = 3;
    private static final int READ_LINE_COUNT = 4;

    private static final char READ_ID_START = '@';
    private static final char READ_ID_BREAK = ' ';
    private static final char READ_ID_DELIM = ':';

    private BufferedWriter createOutputWriter(final String inputFile)
    {
        String fastqFile = inputFile.substring(inputFile.lastIndexOf(File.separator) + 1);
        int extensionIndex = fastqFile.contains("fastq") ? fastqFile.lastIndexOf(".fastq") : fastqFile.lastIndexOf(".fq");
        String outputFile = mOutputDir + fastqFile.substring(0, extensionIndex) + '.' + mOutputId + fastqFile.substring(extensionIndex);

        try
        {
            return outputFile.endsWith(".gz") ? createGzipBufferedWriter(outputFile) : createBufferedWriter(outputFile);
        }
        catch(IOException e)
        {
            MD_LOGGER.error("error creating fastq output file({}): {}", outputFile, e.toString());
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
                        MD_LOGGER.error("invalid entries at line({})", lineCount);

                        for(int i = 0; i < r1ReadBuffer.length; ++i)
                        {
                            MD_LOGGER.error("R1 item {}: {}", i, r1ReadBuffer[i]);
                            MD_LOGGER.error("R2 item {}: {}", i, r2ReadBuffer[i]);
                        }

                        System.exit(1);
                    }

                    readLineCount = 0;
                }

                ++lineCount;

                if(lineCount > 0 && (lineCount % LINE_LOG_COUNT) == 0)
                {
                    MD_LOGGER.info("processed {} lines", lineCount);
                }
            }

            mWriterR1.close();
            mWriterR2.close();
        }
        catch(IOException e)
        {
            MD_LOGGER.error("error reading fastq({}): {}", mFastqFiles, e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }

    private boolean processReadBases(final String[] r1ReadBuffer, final String[] r2ReadBuffer)
    {
        /*
        @A00121:853:H5JJNDSX7:1:1101:1443:1047 1:N:0:GGCACAACCT+CAGGAGTCTA
        GNGAGATGGAGAATTTTCTGGAGATGTCTGAGGAATTTTTTCCTCAGTCTTAAGAGTAAGGTAGGATTGGGCCAGGCATGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCTGAGGCGGACGGATCATGAGGTCAGGAGATCGAG
        +
        F#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
        */

        // data validation
        if(r1ReadBuffer[READ_ITEM_ID].charAt(0) != READ_ID_START || r2ReadBuffer[READ_ITEM_ID].charAt(0) != READ_ID_START)
            return false;

        if(r1ReadBuffer[READ_ITEM_BASES].length() <= mUmiLength || r1ReadBuffer[READ_ITEM_QUALS].length() <= mUmiLength)
            return false;

        int delimIndex = r1ReadBuffer[READ_ITEM_ID].indexOf(READ_ID_BREAK);

        if(delimIndex < 1 || r2ReadBuffer[READ_ITEM_ID].charAt(delimIndex) != READ_ID_BREAK)
            return false;

        String readId1 = r1ReadBuffer[READ_ITEM_ID].substring(1, delimIndex);
        String readId2 = r2ReadBuffer[READ_ITEM_ID].substring(1, delimIndex);

        if(!readId1.equals(readId2))
            return false;

        String umiBases1 = r1ReadBuffer[READ_ITEM_BASES].substring(0, mUmiLength);
        String umiBases2 = r2ReadBuffer[READ_ITEM_BASES].substring(0, mUmiLength);

        // append UMIs to read Id and remove from bases and quals
        String duplexUmiId = umiBases1 + DEFAULT_DUPLEX_UMI_DELIM + umiBases2;
        String newReadId = readId1 + READ_ID_DELIM + duplexUmiId;
        r1ReadBuffer[READ_ITEM_ID] = READ_ID_START + newReadId + r1ReadBuffer[READ_ITEM_ID].substring(delimIndex);
        r2ReadBuffer[READ_ITEM_ID] = READ_ID_START + newReadId + r2ReadBuffer[READ_ITEM_ID].substring(delimIndex);

        r1ReadBuffer[READ_ITEM_BASES] = r1ReadBuffer[READ_ITEM_BASES].substring(mUmiLength);
        r2ReadBuffer[READ_ITEM_BASES] = r2ReadBuffer[READ_ITEM_BASES].substring(mUmiLength);
        r1ReadBuffer[READ_ITEM_QUALS] = r1ReadBuffer[READ_ITEM_QUALS].substring(mUmiLength);
        r2ReadBuffer[READ_ITEM_QUALS] = r2ReadBuffer[READ_ITEM_QUALS].substring(mUmiLength);

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
            MD_LOGGER.error("failed to write output file: {}", e.toString());
            return false;
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(FASTQ_FILES, true, "Fastq file-pair path, separated by delim ','");
        configBuilder.addConfigItem(ConfigItemType.INTEGER, UMI_LENGTH, true, "UMI length", null);
        addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        FastqUmiExtracter fastqUmiExtracter = new FastqUmiExtracter(configBuilder);
        fastqUmiExtracter.run();
    }
}
