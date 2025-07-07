package com.hartwig.hmftools.fastqtools;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASES;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.fastqtools.FastqCommon.FQ_LOGGER;

import org.jetbrains.annotations.NotNull;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Map;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class FastqAnalyser
{
    private final String mFastqFile;
    private final String mOutputFile;
    private final RepeatFrequencies mRepeatFrequencies;

    // repeat sequences can span lines so need to keep the current sequence state
    private StringBuilder mCurrentSequence;
    private char mLastBase;

    // config
    private static final String FASTQ_FILE = "fastq_file";
    private static final String OUTPUT_FILE = "output_file";

    private static final String FASTQ_LINE_START = ">";
    private static final char UNMAPPED_BASE = 'N';
    private static final char NO_BASE = 0;
    private static final int LINE_LOG_COUNT = 1000000;

    public FastqAnalyser(final ConfigBuilder configBuilder)
    {
        mFastqFile = configBuilder.getValue(FASTQ_FILE);
        mOutputFile = configBuilder.getValue(OUTPUT_FILE);
        mRepeatFrequencies = new RepeatFrequencies();
        mCurrentSequence = null;
        mLastBase = NO_BASE;
    }

    public void run()
    {
        if(mOutputFile == null || mFastqFile == null)
            System.exit(1);

        FQ_LOGGER.info("Starting Fastq Analyser with file: {}", mFastqFile);


        long startTimeMs = System.currentTimeMillis();

        analyseFastqFile();

        FQ_LOGGER.info("Fastq analysis complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private void analyseFastqFile()
    {
        try
        {
            BufferedReader fileReader = createBufferedReader(mFastqFile);

            String line = null;
            int lineCount = 0;

            while((line = fileReader.readLine()) != null)
            {
                if(line.startsWith(FASTQ_LINE_START))
                {
                    if(mCurrentSequence != null)
                    {
                        mRepeatFrequencies.registerSequence(mCurrentSequence.toString());
                        mLastBase = NO_BASE;
                        mCurrentSequence = null;
                    }

                    continue;
                }

                processReadBases(line);
                mLastBase = line.charAt(line.length() - 1);

                ++lineCount;

                if(lineCount > 0 && (lineCount % LINE_LOG_COUNT) == 0)
                {
                    FQ_LOGGER.info("processed {} lines", lineCount);
                }
            }

            if(mCurrentSequence != null)
                mRepeatFrequencies.registerSequence(mCurrentSequence.toString());

            // write results
            writeResults();
        }
        catch(IOException e)
        {
            FQ_LOGGER.error("error reading fastq({}): {}", mFastqFile, e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }

    private void processReadBases(final String readBases)
    {
        for(int i = 0; i < readBases.length(); ++i)
        {
            char base = readBases.charAt(i);

            mRepeatFrequencies.registerBase(base);

            // continue building the same sequence
            if(mCurrentSequence != null && mLastBase == base)
            {
                mCurrentSequence.append(base);
                continue;
            }

            // either no current sequence or the repeat has ended
            if(mCurrentSequence != null)
            {
                mRepeatFrequencies.registerSequence(mCurrentSequence.toString());
                mCurrentSequence = null;
            }

            // check for a new repeat starting
            if(base != UNMAPPED_BASE && i < readBases.length() - 1 && readBases.charAt(i + 1) == base)
            {
                mCurrentSequence = new StringBuilder();
                mCurrentSequence.append(base);
                mLastBase = base;
            }
        }
    }

    private void writeResults()
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFile, false);
            writer.write("Sequence,Count");
            writer.newLine();

            for(int i = 0; i < mRepeatFrequencies.baseCounts().length; ++i)
            {
                writer.write(format("%c,%d", DNA_BASES[i], mRepeatFrequencies.baseCounts()[i]));
                writer.newLine();
            }

            for(Map.Entry<String,Integer> entry : mRepeatFrequencies.repeatFrequencies().entrySet())
            {
                String sequence = entry.getKey();
                String seqCode = format("%c%d", sequence.charAt(0), sequence.length());
                writer.write(format("%s,%d", seqCode, entry.getValue()));
                writer.newLine();
            }

            writer.close();
        }
        catch (IOException e)
        {
            FQ_LOGGER.error("failed to create output file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
        configBuilder.addPath(FASTQ_FILE, true, "Fastq file path");
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output fie");

        configBuilder.checkAndParseCommandLine(args);
        setLogLevel(configBuilder);

        FastqAnalyser fastqAnalyser = new FastqAnalyser(configBuilder);
        fastqAnalyser.run();
    }
}
