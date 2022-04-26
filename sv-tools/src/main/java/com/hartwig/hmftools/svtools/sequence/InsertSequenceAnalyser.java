package com.hartwig.hmftools.svtools.sequence;

import static java.lang.Math.round;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.utils.ConfigUtils;

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
    private final Map<String,List<InsertSeqData>> mSampleInsertSeqData;
    private BufferedWriter mWriter;

    private int mMinSearchLength;
    private double mRequiredMatchPerc;

    private List<String> mMatchSequences;

    private static final String INS_SEQ_FILE = "insert_seq_file";
    private static final String SEARCH_LENGTH_MIN = "search_length_min";
    private static final String REQD_MATCH_PERC = "req_match_perc";

    private static final Logger LOGGER = LogManager.getLogger(InsertSequenceAnalyser.class);

    public InsertSequenceAnalyser(final String inputFile, final String outputDir, int minSearch, double reqMatchCount)
    {
        mInputFile = inputFile;
        mOutputDir = outputDir;
        mMinSearchLength = minSearch;
        mRequiredMatchPerc = reqMatchCount;

        mSampleInsertSeqData = Maps.newHashMap();
        mMatchSequences = Lists.newArrayList();

        mWriter = initialiseWriter();
    }

    public void run()
    {
        loadSampleData();

        if(mWriter == null)
            return;

        analyseSampleSequences();
        closeBufferedWriter(mWriter);
    }

    private void analyseSampleSequences()
    {
        int sampleCount = 0;
        for(Map.Entry<String,List<InsertSeqData>> entry : mSampleInsertSeqData.entrySet())
        {
            String sampleId = entry.getKey();
            List<InsertSeqData> insertSeqDataList = entry.getValue();

            LOGGER.trace("sample({}) checking {} SGLs", sampleId, insertSeqDataList.size());

            int index = 0;
            while(index < insertSeqDataList.size())
            {
                InsertSeqData first = insertSeqDataList.get(index);

                boolean matched = false;

                for(int nextIndex = index + 1; nextIndex < insertSeqDataList.size(); ++nextIndex)
                {
                    InsertSeqData next = insertSeqDataList.get(nextIndex);

                    MatchedBaseData matchedBaseData = getMatchedBaseData(first, next);

                    if(matchedBaseData != null)
                    {
                        writeMatchData(sampleId, first, next, matchedBaseData);
                        matched = true;
                        insertSeqDataList.remove(nextIndex);
                        break;
                    }
                }

                if(matched)
                    insertSeqDataList.remove(index);
                else
                    ++index;
            }

            ++sampleCount;

            if((sampleCount % 100) == 0)
            {
                LOGGER.info("processed {} samples", sampleCount);
            }
        }

        LOGGER.info("insert sequence analysis complete");
    }

    private MatchedBaseData getMatchedBaseData(final InsertSeqData sgl1, final InsertSeqData sgl2)
    {
        // scenarios:
        // different orientations: the -ve orientation will read the insert sequence as reverse compliment and in reverse order
        // same orientations: the second SGL will read the insert sequence as forward compliment but in reverse order

        String firstSequence = sgl1.InsertSeq;
        String secondSequence = sgl2.InsertSeq;

        if(sgl1.Orientation == sgl2.Orientation)
        {
            secondSequence = Nucleotides.reverseStrandBases(secondSequence);
        }
        else
        {
            //
            // secondSequence = reverseString(secondSequence);
        }

        MatchedBaseData matchData = findMatchData(firstSequence, secondSequence);

        if(matchData != null)
            return matchData;

        matchData = findMatchData(secondSequence, firstSequence);

        return matchData;
    }

    private static final int START_SEARCH_LENGTH = 20;

    private MatchedBaseData findMatchData(final String seq1, final String seq2)
    {
        int matchStart1 = -1;
        int matchStart2 = -1;

        for(int i = 0; i < seq2.length() - START_SEARCH_LENGTH; ++i)
        {
            String searchStr = seq2.substring(i, i + START_SEARCH_LENGTH);

            int matchIndex = seq1.indexOf(searchStr);

            if(matchIndex >= 0)
            {
                matchStart1 = matchIndex;
                matchStart2 = i;
                break;
            }
        }

        if(matchStart1 < 0)
            return null;

        int matchedBases = START_SEARCH_LENGTH;
        int misMatches = 0;
        int currentMismatches = 0;

        for(int index2 = matchStart2 + matchedBases; index2 < seq2.length(); ++index2)
        {
            int index1 = matchStart1 + matchedBases;
            if(index1 >= seq1.length())
                break;

            if(seq1.charAt(index1) == seq2.charAt(index2))
            {
                ++matchedBases;
                currentMismatches = 0;
                continue;
            }

            ++misMatches;
            ++currentMismatches;

            if(currentMismatches >= 2)
            {
                misMatches -= currentMismatches;
                break;
            }

            if(misMatches / (double)(matchedBases + misMatches) > (1 - mRequiredMatchPerc))
            {
                --misMatches;
                break;
            }
        }

        return matchedBases >= mMinSearchLength ? new MatchedBaseData(matchedBases + misMatches, matchedBases) : null;
    }

    private StructuralVariantType getSvType(final InsertSeqData sgl1, final InsertSeqData sgl2)
    {
        if(!sgl1.Chromosome.equals(sgl2.Chromosome))
            return BND;

        if(sgl1.Orientation == sgl2.Orientation)
        {
            return INV;
        }

        if(sgl1.Position < sgl2.Position)
            return sgl1.Orientation == POS_ORIENT ? DEL : DUP;
        else
            return sgl2.Orientation == POS_ORIENT ? DEL : DUP;
    }

    private BufferedWriter initialiseWriter()
    {
        if(mOutputDir == null)
            return null;

        try
        {
            final String outputFileName = mOutputDir + "insert_seq_matches.csv";

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("SampleId,SvId1,SvId2,VcfId1,VcdId2,Chr1,Pos1,Orient1,Chr2,Pos2,Orient2");
            writer.write(",SvType,OverlapBases,MatchedBases,SeqLen1,SeqLen2");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            LOGGER.error("Failed to write insert sequence output CSV file: {}", e.toString());
            return null;
        }
    }

    private void writeMatchData(final String sampleId, final InsertSeqData first, final InsertSeqData next, final MatchedBaseData matchedBases)
    {
        try
        {
            StructuralVariantType svType = getSvType(first, next);

            mWriter.write(String.format("%s,%d,%d,%s,%s,%s,%d,%d,%s,%d,%d",
                    sampleId, first.SvId, next.SvId, first.VcfId, next.VcfId,
                    first.Chromosome, first.Position, first.Orientation,
                    next.Chromosome, next.Position, next.Orientation));

            mWriter.write(String.format(",%s,%d,%d,%d,%d",
                    svType, matchedBases.OverlapBases, matchedBases.MatchedBases, first.InsertSeq.length(), next.InsertSeq.length()));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            LOGGER.error("Failed to write insert sequence matches: {}", e.toString());
        }
    }

    private void loadSampleData()
    {
        if(mInputFile.isEmpty())
            return;

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(mInputFile));
            lines.remove(0);

            List<InsertSeqData> sampleDataList = null;
            String currentSample = "";

            for(String line : lines)
            {
                InsertSeqData insData = InsertSeqData.fromCsv(line);

                if(insData.InsertSeq.length() < mMinSearchLength)
                    continue;

                if(!currentSample.equals(insData.SampleId))
                {
                    sampleDataList = Lists.newArrayList();
                    currentSample = insData.SampleId;
                    mSampleInsertSeqData.put(insData.SampleId, sampleDataList);
                }

                sampleDataList.add(insData);
            }

            LOGGER.info("loaded samples({}) total insert sequences from CSV file({})",
                    mSampleInsertSeqData.size(), mSampleInsertSeqData.values().stream().mapToInt(x -> x.size()).sum(), mInputFile);

        }
        catch(IOException e)
        {
            LOGGER.error("failed to read insert sequence CSV file({}): {}", mInputFile, e.toString());
        }
    }

    private class MatchedBaseData
    {
        public int OverlapBases;
        public int MatchedBases;

        public MatchedBaseData(final int overlapBases, final int matchedBases)
        {
            OverlapBases = overlapBases;
            MatchedBases = matchedBases;
        }
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

    public static InsertSequenceAnalyser from (final CommandLine cmd)
    {
        return new InsertSequenceAnalyser(
                cmd.getOptionValue(INS_SEQ_FILE),
                parseOutputDir(cmd),
                Integer.parseInt(cmd.getOptionValue(SEARCH_LENGTH_MIN, "50")),
                Double.parseDouble(cmd.getOptionValue(REQD_MATCH_PERC, "0.95")));
    }

    public static void addCmdLineArgs(final Options options)
    {
        options.addOption(INS_SEQ_FILE, true, "File with sample insert sequence data");
        options.addOption(SEARCH_LENGTH_MIN, true, "Sequence min length");
        options.addOption(REQD_MATCH_PERC, true, "Number of instances of search string to record it");
        addOutputDir(options);
        addLoggingOptions(options);
    }

    /*

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

                    if(matches.size() >= mRequiredMatchPerc)
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
    */

}
