package com.hartwig.hmftools.fastqtools.umi;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadDelimitedIdFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createGzipBufferedWriter;
import static com.hartwig.hmftools.fastqtools.FastqCommon.APP_NAME;
import static com.hartwig.hmftools.fastqtools.FastqCommon.FQ_LOGGER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class FastqUmiExtracter
{
    private final UmiConfig mConfig;

    private BufferedWriter mWriterR1;
    private BufferedWriter mWriterR2;

    private final List<String> mKnownUmis;
    private final int mMaxKnownUmiLength;

    private static final String FASTQ_FILES_DELIM = ";";
    private static final int LINE_LOG_COUNT = 1000000;

    public FastqUmiExtracter(final ConfigBuilder configBuilder)
    {
        mConfig = new UmiConfig(configBuilder);
        mKnownUmis = Lists.newArrayList();

        if(mConfig.KnownUmiFile != null)
            loadKnownUmis();

        mMaxKnownUmiLength = mKnownUmis.stream().mapToInt(x -> x.length()).max().orElse(0);
    }

    private void loadKnownUmis()
    {
        List<String> knownUmis = loadDelimitedIdFile(mConfig.KnownUmiFile, "KnownUmi", CSV_DELIM);
        knownUmis.forEach(x -> mKnownUmis.add(x));

        Collections.sort(mKnownUmis, Comparator.comparingInt(x -> -x.length()));

        if(!knownUmis.isEmpty())
        {
            FQ_LOGGER.info("loaded {} known UMIs from {}", mKnownUmis.size(), mConfig.KnownUmiFile);
        }
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

        FQ_LOGGER.info("extraction complete, mins({})", runTimeMinsStr(startTimeMs));
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
        /*
        @A00121:853:H5JJNDSX7:1:1101:1443:1047 1:N:0:GGCACAACCT+CAGGAGTCTA
        GNGAGATGGAGAATTTTCTGGAGATGTCTGAGGAATTTTTTCCTCAGTCTTAAGAGTAAGGTAGGATTGGGCCAGGCATGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCTGAGGCGGACGGATCATGAGGTCAGGAGATCGAG
        +
        F#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
        */

        // data validation
        if(r1ReadBuffer[READ_ITEM_ID].charAt(0) != READ_ID_START || r2ReadBuffer[READ_ITEM_ID].charAt(0) != READ_ID_START)
            return false;
        
        int minReadLength = max(mConfig.AdapterUmiLength, mConfig.UmiLength);

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
            adjustWithAdapter(readId1, readId2, delimIndex, r1ReadBuffer,  r2ReadBuffer);
        }
        else if(!mKnownUmis.isEmpty())
        {
            adjustWithKnownUmi(readId1, delimIndex, r1ReadBuffer, r2ReadBuffer);
        }
        else
        {
            adjustWithFixedUmi(readId1, delimIndex, r1ReadBuffer, r2ReadBuffer);
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

    private void adjustWithKnownUmi(
            final String readId1, final int delimIndex, final String[] r1ReadBuffer, final String[] r2ReadBuffer)
    {
        // find the longest matching UMI

        String umiBases1 = findKnownUmiMatch(r1ReadBuffer[READ_ITEM_BASES], false);
        String umiBases2 = findKnownUmiMatch(r2ReadBuffer[READ_ITEM_BASES], true);

        // make UMIs of the form 4+5, where integers are the length of the UMIs
        int umiLength1 = umiBases1.length();
        int umiLength2 = umiBases2.length();

        // append UMIs to read Id and remove from bases and quals
        String duplexUmiId = umiLength1 + mConfig.UmiDelim + umiLength2;
        String newReadId = readId1 + READ_ID_DELIM + duplexUmiId;
        r1ReadBuffer[READ_ITEM_ID] = READ_ID_START + newReadId + r1ReadBuffer[READ_ITEM_ID].substring(delimIndex);
        r2ReadBuffer[READ_ITEM_ID] = READ_ID_START + newReadId + r2ReadBuffer[READ_ITEM_ID].substring(delimIndex);

        r1ReadBuffer[READ_ITEM_BASES] = r1ReadBuffer[READ_ITEM_BASES].substring(umiLength1);
        r2ReadBuffer[READ_ITEM_BASES] = r2ReadBuffer[READ_ITEM_BASES].substring(umiLength2);
        r1ReadBuffer[READ_ITEM_QUALS] = r1ReadBuffer[READ_ITEM_QUALS].substring(umiLength1);
        r2ReadBuffer[READ_ITEM_QUALS] = r2ReadBuffer[READ_ITEM_QUALS].substring(umiLength2);
    }

    private void adjustWithFixedUmi(
            final String readId1, final int delimIndex, final String[] r1ReadBuffer, final String[] r2ReadBuffer)
    {
        String umiBases1 = r1ReadBuffer[READ_ITEM_BASES].substring(0, mConfig.UmiLength);
        String umiBases2 = r2ReadBuffer[READ_ITEM_BASES].substring(0, mConfig.UmiLength);

        // append UMIs to read Id and remove from bases and quals
        String duplexUmiId = umiBases1 + mConfig.UmiDelim + umiBases2;
        String newReadId = readId1 + READ_ID_DELIM + duplexUmiId;
        r1ReadBuffer[READ_ITEM_ID] = READ_ID_START + newReadId + r1ReadBuffer[READ_ITEM_ID].substring(delimIndex);
        r2ReadBuffer[READ_ITEM_ID] = READ_ID_START + newReadId + r2ReadBuffer[READ_ITEM_ID].substring(delimIndex);

        r1ReadBuffer[READ_ITEM_BASES] = r1ReadBuffer[READ_ITEM_BASES].substring(mConfig.UmiLength);
        r2ReadBuffer[READ_ITEM_BASES] = r2ReadBuffer[READ_ITEM_BASES].substring(mConfig.UmiLength);
        r1ReadBuffer[READ_ITEM_QUALS] = r1ReadBuffer[READ_ITEM_QUALS].substring(mConfig.UmiLength);
        r2ReadBuffer[READ_ITEM_QUALS] = r2ReadBuffer[READ_ITEM_QUALS].substring(mConfig.UmiLength);
    }

    private void adjustWithAdapter(
            final String readId1, final String readId2, final int delimIndex,
            final String[] r1ReadBuffer, final String[] r2ReadBuffer)
    {
        String adapterUmiBases1 = r1ReadBuffer[READ_ITEM_BASES].substring(0, mConfig.AdapterUmiLength);

        String umiBases1 = adapterUmiBases1.substring(0, mConfig.UmiLength);

        // append UMIs to read Id and remove from bases and quals
        String newReadId = readId1 + READ_ID_DELIM + umiBases1;
        r1ReadBuffer[READ_ITEM_ID] = READ_ID_START + newReadId + r1ReadBuffer[READ_ITEM_ID].substring(delimIndex);
        r2ReadBuffer[READ_ITEM_ID] = READ_ID_START + newReadId + r2ReadBuffer[READ_ITEM_ID].substring(delimIndex);

        r1ReadBuffer[READ_ITEM_BASES] = r1ReadBuffer[READ_ITEM_BASES].substring(mConfig.AdapterUmiLength);
        r1ReadBuffer[READ_ITEM_QUALS] = r1ReadBuffer[READ_ITEM_QUALS].substring(mConfig.AdapterUmiLength);

        // the R2 read may have the reversed adapter+UMI sequence at the end
        String adapterUmiBases2Start = r2ReadBuffer[READ_ITEM_BASES].substring(0, mConfig.AdapterUmiLength);
        int read2Length = r2ReadBuffer[READ_ITEM_BASES].length();
        String adapterUmiBases2End = r2ReadBuffer[READ_ITEM_BASES].substring(read2Length - mConfig.AdapterUmiLength);

        int adapterSeqIndex = adapterUmiBases2Start.indexOf(mConfig.AdapterSequence);

        if(adapterSeqIndex >= 0)
        {
            // trim from start
            int trimIndex = adapterSeqIndex + mConfig.AdapterSequence.length();
            r2ReadBuffer[READ_ITEM_BASES] = r2ReadBuffer[READ_ITEM_BASES].substring(trimIndex);
            r2ReadBuffer[READ_ITEM_QUALS] = r2ReadBuffer[READ_ITEM_QUALS].substring(trimIndex);
        }
        else
        {
            adapterSeqIndex = adapterUmiBases2End.indexOf(mConfig.AdapterSequenceReversed);

            if(adapterSeqIndex >= 0)
            {
                int trimIndex = read2Length - mConfig.AdapterUmiLength + adapterSeqIndex;
                r2ReadBuffer[READ_ITEM_BASES] = r2ReadBuffer[READ_ITEM_BASES].substring(0, trimIndex);
                r2ReadBuffer[READ_ITEM_QUALS] = r2ReadBuffer[READ_ITEM_QUALS].substring(0, trimIndex);
            }
        }
    }

    private String findKnownUmiMatch(final String readBases, boolean requiresReverse)
    {
        for(String umi : mKnownUmis)
        {
            String readUmi = readBases.substring(0, umi.length());

            if(requiresReverse)
                readUmi = Nucleotides.reverseComplementBases(readUmi);

            if(readUmi.equals(umi))
                return umi;
        }

        // check for a 1-base mismatch
        for(String umi : mKnownUmis)
        {
            String readUmi = readBases.substring(0, umi.length());

            if(requiresReverse)
                readUmi = Nucleotides.reverseComplementBases(readUmi);

            if(!exceedsUmiDiff(readUmi, umi, mConfig.KnownUmiBaseDiff))
                return umi;
        }

        return "";
    }

    private static boolean exceedsUmiDiff(final String first, final String second, int permittedDiff)
    {
        if(first.length() != second.length())
            return true;

        short diffs = 0;
        for(short i = 0; i < first.length(); ++i)
        {
            if(first.charAt(i) != second.charAt(i))
            {
                ++diffs;

                if(diffs > permittedDiff)
                    return true;
            }
        }

        return false;
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
