package com.hartwig.hmftools.svtools.sequence;

import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class InsertSequenceAnalyser
{
    private final String mInputFile;
    private final String mOutputDir;
    private List<InsertSeqData> mInsertSeqData;
    private Map<String,List<InsertSeqData>> mSampleInsertSeqData;
    private BufferedWriter mWriter;

    private int mMinSearchLength;
    private int mMaxSearchLength;
    private int mRequiredMatchCount;

    private List<String> mMatchSequences;

    private static final String INS_SEQ_INPUT_FILE = "ins_seq_input_file";
    private static final String SEARCH_LENGTH_MIN = "search_length_min";
    private static final String SEARCH_LENGTH_MAX = "search_length_max";
    private static final String REQD_MATCH_COUNT = "req_match_count";

    private static final Logger LOGGER = LogManager.getLogger(InsertSequenceAnalyser.class);

    public InsertSequenceAnalyser(final String inputFile, final String outputDir, int minSearch, int maxSearch, int reqMatchCount)
    {
        mInputFile = inputFile;
        mOutputDir = outputDir;
        mMinSearchLength = minSearch;
        mMaxSearchLength = maxSearch;
        mRequiredMatchCount = reqMatchCount;

        mInsertSeqData = Lists.newArrayList();
        mSampleInsertSeqData = Maps.newHashMap();
        mMatchSequences = Lists.newArrayList();

        mWriter = null;
    }

    public static InsertSequenceAnalyser from (final CommandLine cmd)
    {
        return new InsertSequenceAnalyser(
                cmd.getOptionValue(INS_SEQ_INPUT_FILE),
                parseOutputDir(cmd),
                Integer.parseInt(cmd.getOptionValue(SEARCH_LENGTH_MAX, "-1")),
                Integer.parseInt(cmd.getOptionValue(SEARCH_LENGTH_MIN, "1")),
                Integer.parseInt(cmd.getOptionValue(REQD_MATCH_COUNT, "10")));
    }

    public static void addCmdLineArgs(final Options options)
    {
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(INS_SEQ_INPUT_FILE, true, "File with sample insert sequence data");
        options.addOption(SEARCH_LENGTH_MAX, true, "Sequence max length");
        options.addOption(SEARCH_LENGTH_MIN, true, "Sequence min length");
        options.addOption(REQD_MATCH_COUNT, true, "Number of instances of search string to record it");
    }

    public final List<String> getMatchSequences() { return mMatchSequences; }

    @VisibleForTesting
    public void addSequenceData(final List<InsertSeqData> seqDataList)
    {
        mInsertSeqData.addAll(seqDataList);

        List<InsertSeqData> sampleDataList = null;
        String currentSample = "";

        for(final InsertSeqData isData : seqDataList)
        {
            if(!currentSample.equals(isData.SampleId))
            {
                sampleDataList = Lists.newArrayList();
                currentSample = isData.SampleId;
                mSampleInsertSeqData.put(isData.SampleId, sampleDataList);
            }

            sampleDataList.add(isData);
        }
    }

    public void run()
    {
        loadSampleData();
        analyseCommonSequences();
        analyseSampleSequences();
        closeBufferedWriter(mWriter);
    }

    private static final int COL_SAMPLE_ID = 0;
    private static final int COL_SV_ID = 1;
    private static final int COL_INS_SEQ = 2;

    private void loadSampleData()
    {
        if(mInputFile.isEmpty())
            return;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(mInputFile));

            String line = fileReader.readLine(); // skip header

            List<InsertSeqData> sampleDataList = null;
            String currentSample = "";

            while ((line = fileReader.readLine()) != null)
            {
                // parse CSV data
                String[] items = line.split(",", -1);

                if(items.length < COL_INS_SEQ+1)
                {
                    LOGGER.warn("invalid input string: {}", line);
                    continue;
                }

                final String sampleId = items[COL_SAMPLE_ID];
                int svId = Integer.parseInt(items[COL_SV_ID]);
                final String insertSeq = items[COL_INS_SEQ];

                InsertSeqData isData = new InsertSeqData(sampleId, svId, insertSeq);
                mInsertSeqData.add(isData);

                if(!currentSample.equals(sampleId))
                {
                    sampleDataList = Lists.newArrayList();
                    currentSample = sampleId;
                    mSampleInsertSeqData.put(sampleId, sampleDataList);
                }

                sampleDataList.add(isData);
            }

            LOGGER.info("loaded {} sample insert sequences from CSV file({})", mInsertSeqData.size(), mInputFile);

        }
        catch(IOException e)
        {
            LOGGER.error("failed to read insert sequence CSV file({}): {}", mInputFile, e.toString());
        }
    }

    private void analyseSampleSequences()
    {

    }

    private void analyseCommonSequences()
    {
        for(int i = 0; i < mInsertSeqData.size() - 1; ++i)
        {
            if(i > 0 && (i % 100) == 0)
            {
                LOGGER.info("processed {} sequences", i);
            }

            final InsertSeqData isData = mInsertSeqData.get(i);
            final String insertSeq = isData.InsertSeq;
            int isLength = insertSeq.length();
            int maxLength = mMaxSearchLength > 0 ? mMaxSearchLength : isLength;

            // search for the longest match first
            for(int len = maxLength; len >= mMinSearchLength; --len)
            {
                for(int strIndex = 0; strIndex < isLength; ++strIndex)
                {
                    if(strIndex + len > isLength)
                        break;

                    String searchStr = insertSeq.substring(strIndex, strIndex+len);

                    if(searchStr.isEmpty())
                        continue;

                    if(hasMatchedSequence(searchStr))
                        continue;

                    List<InsertSeqData> matches = findSequenceMatches(searchStr, mInsertSeqData, i+1);

                    if(matches.isEmpty())
                        continue;

                    if(matches.size() >= mRequiredMatchCount)
                    {
                        mMatchSequences.add(searchStr);
                        logMatchSequence(false, isData, matches, searchStr);

                        // move forward to the start of the next sequence by the min required overlap amount
                        strIndex += (int)round(len * (1 - MIN_MATCH_PERCENT));
                    }
                }
            }
        }
    }

    private static double MIN_MATCH_PERCENT = 0.90;

    private boolean hasMatchedSequence(final String searchStr)
    {
        if(mMatchSequences.isEmpty())
            return false;

        if(mMatchSequences.contains(searchStr))
            return true;

        // check for high levels of overlap
        for(final String str : mMatchSequences)
        {
            final String shortStr = searchStr.length() < str.length() ? searchStr : str;
            final String longStr = searchStr.length() >= str.length() ? searchStr : str;
            int shortLen = shortStr.length();

            int minLength = (int) round(MIN_MATCH_PERCENT * shortLen);

            for (int strIndex = 0; strIndex < shortLen; ++strIndex)
            {
                for (int len = minLength; len <= shortLen; ++len)
                {
                    if (strIndex + len > shortLen)
                        break;

                    final String subStr = shortStr.substring(strIndex, strIndex + len);

                    if (longStr.contains(subStr))
                        return true;
                }
            }
        }

        return false;
    }

    private List<InsertSeqData> findSequenceMatches(final String searchStr, List<InsertSeqData> searchList, Integer startIndex)
    {
        final List<InsertSeqData> matches = Lists.newArrayList();

        int i = startIndex != null ? startIndex : 0;

        for(; i < mInsertSeqData.size(); ++i)
        {
            InsertSeqData isData = searchList.get(i);

            if(isData.InsertSeq.contains(searchStr))
                matches.add(isData);
        }

        return matches;
    }

    private final String reverseString(final String str)
    {
        String reverse = "";

        for(int i = str.length() - 1; i >= 0; --i)
        {
            reverse += str.charAt(i);
        }

        return reverse;
    }

    private void logMatchSequence(boolean specificSample, final InsertSeqData isData, final List<InsertSeqData> matches, final String matchSequence)
    {
        if(mOutputDir.isEmpty())
            return;

        try
        {
            if(mWriter == null)
            {
                final String outputFileName = mOutputDir + "LNX_INS_SEQUENCES.csv";

                mWriter = createBufferedWriter(outputFileName, false);
                mWriter.write("SampleId,SampleCount,MatchCount,InsertSequence");
                mWriter.newLine();
            }

            int sampleCount = 1;
            String sampleId = isData.SampleId;

            if(!specificSample)
            {
                // determine unique sample count
                List<String> sampleIds = Lists.newArrayList(isData.SampleId);

                for (final InsertSeqData otherData : matches)
                {
                    if (!sampleIds.contains(otherData.SampleId))
                        sampleIds.add(otherData.SampleId);
                }

                sampleCount = sampleIds.size();
                sampleId = "ALL";
            }

            // +1 is to account for the data which initiates the search
            mWriter.write(String.format("%s,%d,%d,%s", sampleId, sampleCount, matches.size() + 1, matchSequence));
            mWriter.newLine();
        }
        catch(IOException e)
        {
            LOGGER.error("Failed to write insert sequence output CSV file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        addCmdLineArgs(options);

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        Configurator.setRootLevel(Level.DEBUG);

        InsertSequenceAnalyser insSeqAnalyser = InsertSequenceAnalyser.from(cmd);
        insSeqAnalyser.run();
    }

}
