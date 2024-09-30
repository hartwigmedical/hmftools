package com.hartwig.hmftools.bamtools.fastqbimodalcollapse;

import static java.lang.Math.min;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.codon.Nucleotides.baseIndex;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.fastq.FastqRecord;

// TODO: Bimodal -> Biomodal
// TODO: Multithread
// TODO: Break out into separate classes
// TODO: Reduce reimplementation of common classes.
// TODO: Cleanup
public class FastqBimodalCollapse
{
    private static final char MISSING_BASE = 'N';
    private static final char MISMATCH_BASE = 'X';
    private static final char MODC_BASE = 'c';

    private static final int ROLLING_WINDOW_SIZE = 20;
    private static final int MIN_GOOD_ROLLING_COUNT = 18;
    private static final int NO_HAIRPIN = Integer.MAX_VALUE;

    private static final float HAIRPIN_MATCH_TARGET = 0.9f;
    private static final int MIN_HAIRPIN_MATCH_LENGTH = 3;
    private static final String FORWARD_HAIRPIN = "AATGACGATGCGTTCGAGCATCGTTATT";
    private static final String REVERSE_HAIRPIN = "AATACCGATGCTCGAACGCATCGTCATT";

    private final FastqBimodalCollapseConfig mConfig;

    private int mFastqPairCount;

    public FastqBimodalCollapse(final FastqBimodalCollapseConfig config)
    {
        mConfig = config;

        mFastqPairCount = 0;
    }

    @Nullable
    public FastqRecord nextFastqRecord(final BufferedReader reader) throws IOException
    {
        String rawReadName = reader.readLine();
        if(rawReadName == null)
        {
            return null;
        }

        // TODO: Handle tags.
        String readName = rawReadName.split("\\s+", 2)[0].substring(1);

        String readString = reader.readLine();
        String qualityHeader = reader.readLine();
        String baseQualityString = reader.readLine();
        if(baseQualityString == null)
        {
            throw new RuntimeException("Partial fastq record found.");
        }

        return new FastqRecord(readName, readString, qualityHeader, baseQualityString);
    }

    public void run()
    {
        BT_LOGGER.info("starting FastqBimodalCollapse");
        long startTimeMs = System.currentTimeMillis();

        String statFilePath = Paths.get(mConfig.OutputDir, mConfig.SampleId + ".stats.tsv").toString();

        // TODO: use fastq reader?
        // TODO: read gzipped files.
        try(BufferedReader fastq1Reader = new BufferedReader(new FileReader(mConfig.Fastq1Path));
                BufferedReader fastq2Reader = new BufferedReader(new FileReader(mConfig.Fastq2Path));
                BufferedWriter statWriter = new BufferedWriter(new FileWriter(statFilePath)))
        {
            statWriter.write(Arrays.stream(STAT_HEADERS).collect(Collectors.joining("\t")));
            statWriter.newLine();

            FastqRecord fastq1 = nextFastqRecord(fastq1Reader);
            FastqRecord fastq2 = nextFastqRecord(fastq2Reader);
            while(fastq1 != null && fastq2 != null)
            {
                mFastqPairCount++;

                if(!fastq1.getReadName().equals(fastq2.getReadName()))
                {
                    throw new RuntimeException("Fastq read name mismatch.");
                }

                // TODO: We just trim down?
                //                if(fastq1.getReadLength() != fastq2.getReadLength())
                //                {
                //                    mReadLengthMismatches++;
                //                }

                processFastqPair(statWriter, fastq1, fastq2);

                fastq1 = nextFastqRecord(fastq1Reader);
                fastq2 = nextFastqRecord(fastq2Reader);
            }

            if(fastq1 != null || fastq2 != null)
            {
                throw new RuntimeException("Fastq file terminated early.");
            }
        }
        catch(FileNotFoundException e)
        {
            throw new RuntimeException(e);
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }

        // TODO: Add extra info.
        BT_LOGGER.info("FastqBimodalCollapse complete, mins({}) fastqPairCount({})", runTimeMinsStr(startTimeMs), mFastqPairCount);
    }

    private static char getConsensusBase(char base1, char base2)
    {
        // TODO: Fix logic for N's, but for now this doesn't affect the stats.
        if(baseIndex(base1) == -1 && base1 != MISSING_BASE)
        {
            throw new RuntimeException("Invalid base1");
        }

        if(baseIndex(base2) == -1 && base2 != MISSING_BASE)
        {
            throw new RuntimeException("Invalid base2");
        }

        if(base1 == MISSING_BASE)
        {
            return base2;
        }

        if(base2 == MISSING_BASE)
        {
            return base1;
        }

        if(base2 == 'G')
        {
            return MISMATCH_BASE;
        }

        if(base1 == base2)
        {
            return base1 == 'C' ? MODC_BASE : base1;
        }

        if(base1 == 'G')
        {
            return base2 == 'A' ? 'G' : MISMATCH_BASE;
        }

        if(base1 == 'T')
        {
            return base2 == 'C' ? 'C' : MISMATCH_BASE;
        }

        return MISMATCH_BASE;
    }

    private static String getConsensusRead(final FastqRecord fastq1, final FastqRecord fastq2)
    {
        StringBuilder consensusRead = new StringBuilder();
        int length = min(fastq1.getReadLength(), fastq2.getReadLength());
        for(int i = 0; i < length; i++)
        {
            char base1 = fastq1.getReadString().charAt(i);
            char base2 = fastq2.getReadString().charAt(i);

            consensusRead.append(getConsensusBase(base1, base2));
        }

        return consensusRead.toString();
    }

    public static class HairpinInfo
    {
        public final int StartPos;
        public final int Length;
        public final int MatchCount;

        public HairpinInfo(final int startPos, final int length, final int matchCount)
        {
            StartPos = startPos;
            Length = length;
            MatchCount = matchCount;
        }
    }

    private static HairpinInfo findHairpin(final String read, final String hairpinSequence)
    {
        // full hairpin match
        int minRequiredMatches = (int) Math.round(Math.ceil(hairpinSequence.length() * HAIRPIN_MATCH_TARGET));
        int maxMismatches = hairpinSequence.length() - minRequiredMatches;
        for(int i = 0; i <= read.length() - hairpinSequence.length(); i++)
        {
            int mismatchCount = 0;
            for(int j = 0; j < hairpinSequence.length(); j++)
            {
                char readBase = read.charAt(i + j);
                char hairpinBase = hairpinSequence.charAt(j);
                if(readBase != hairpinBase)
                {
                    mismatchCount++;
                }

                if(mismatchCount > maxMismatches)
                {
                    break;
                }
            }

            if(mismatchCount > maxMismatches)
            {
                continue;
            }

            return new HairpinInfo(i + 1, hairpinSequence.length(), hairpinSequence.length() - mismatchCount);
        }

        // partial hairpin match at the end
        for(int hairpinPrefixLength = min(hairpinSequence.length() - 1, read.length()); hairpinPrefixLength >= MIN_HAIRPIN_MATCH_LENGTH; hairpinPrefixLength--)
        {
            minRequiredMatches = (int) Math.round(Math.ceil(hairpinPrefixLength * HAIRPIN_MATCH_TARGET));
            maxMismatches = hairpinPrefixLength - minRequiredMatches;

            int mismatchCount = 0;
            for(int j = 0; j < hairpinPrefixLength; j++)
            {
                char readBase = read.charAt(read.length() - hairpinPrefixLength + j);
                char hairpinBase = hairpinSequence.charAt(j);
                if(readBase != hairpinBase)
                {
                    mismatchCount++;
                }

                if(mismatchCount > maxMismatches)
                {
                    break;
                }
            }

            if(mismatchCount > maxMismatches)
            {
                continue;
            }

            return new HairpinInfo(read.length() - hairpinPrefixLength + 1, hairpinPrefixLength, hairpinPrefixLength - mismatchCount);
        }

        return new HairpinInfo(-1, 0, 0);
    }

    private static final String STAT_DELIMITER = "\t";
    private static final String[] STAT_HEADERS = {
            "read_name",
            "read1_length",
            "read2_length",
            "read1",
            "qual1",
            "read2",
            "qual2",
            "hairpin1_start_pos",
            "hairpin1_length",
            "hairpin1_match_count",
            "hairpin2_start_pos",
            "hairpin2_length",
            "hairpin2_match_count"
    };

    private static void processFastqPair(final BufferedWriter writer, final FastqRecord fastq1, final FastqRecord fastq2) throws IOException
    {
        HairpinInfo hairpin1 = findHairpin(fastq1.getReadString(), FORWARD_HAIRPIN);
        HairpinInfo hairpin2 = findHairpin(fastq2.getReadString(), REVERSE_HAIRPIN);

        StringJoiner statLine = new StringJoiner(STAT_DELIMITER);
        statLine.add(fastq1.getReadName());
        statLine.add(String.valueOf(fastq1.getReadLength()));
        statLine.add(String.valueOf(fastq2.getReadLength()));
        statLine.add(fastq1.getReadString());
        statLine.add(fastq1.getBaseQualityString());
        statLine.add(fastq2.getReadString());
        statLine.add(fastq2.getBaseQualityString());
        statLine.add(String.valueOf(hairpin1.StartPos));
        statLine.add(String.valueOf(hairpin1.Length));
        statLine.add(String.valueOf(hairpin1.MatchCount));
        statLine.add(String.valueOf(hairpin2.StartPos));
        statLine.add(String.valueOf(hairpin2.Length));
        statLine.add(String.valueOf(hairpin2.MatchCount));

        writer.write(statLine.toString());
        writer.newLine();
    }

//    private static final String STAT_DELIMITER = "\t";
//    private static final String[] STAT_HEADERS = {
//            "read_name",
//            "read1_length",
//            "read2_length",
//            "consensus_count",
//            "missing_count",
//            "mismatch_count",
//            "modc_with_g_count",
//            "modc_without_g_count",
//            "exact_matches",
//            "exact_mismatches",
//            "hairpin_start",
//            "read1",
//            "qual1",
//            "read2",
//            "qual2",
//            "consensus_read"
//    };
//
//    private static void processFastqPair(final BufferedWriter writer, final FastqRecord fastq1, final FastqRecord fastq2) throws IOException
//    {
//        String consensusRead = getConsensusRead(fastq1, fastq2);
//        int hairpinStart = findHairpinStart(fastq1, fastq2, consensusRead);
//
//        int missingCount = 0;
//        int consensusCount = 0;
//        int mismatchCount = 0;
//        int modCWithG = 0;
//        int modCWithoutG = 0;
//        int exactMatches = 0;
//        int exactMismatches = 0;
//        for(int i = 0; i < min(consensusRead.length(), hairpinStart); i++)
//        {
//            char consensusBase = consensusRead.charAt(i);
//            char base1 = fastq1.getReadString().charAt(i);
//            char base2 = fastq2.getReadString().charAt(i);
//
//            if(base1 == MISSING_BASE || base2 == MISSING_BASE)
//            {
//                missingCount++;
//            }
//            else if(baseIndex(consensusBase) != -1)
//            {
//                consensusCount++;
//            }
//            else if(consensusBase == MISMATCH_BASE)
//            {
//                mismatchCount++;
//            }
//            else if(consensusBase == MODC_BASE)
//            {
//                Character nextBase1 = i < fastq1.getReadString().length() - 1 ? fastq1.getReadString().charAt(i + 1) : null;
//                if(nextBase1 != null && nextBase1.equals('G'))
//                {
//                    modCWithG++;
//                }
//                else
//                {
//                    modCWithoutG++;
//                }
//            }
//            else
//            {
//                throw new RuntimeException("Unreachable");
//            }
//
//
//            if(base1 == MISSING_BASE || base2 == MISSING_BASE)
//            {
//                continue;
//            }
//
//            if(base1 == base2)
//            {
//                exactMatches++;
//                continue;
//            }
//
//            exactMismatches++;
//        }
//
//        StringJoiner statLine = new StringJoiner(STAT_DELIMITER);
//        statLine.add(fastq1.getReadName());
//        statLine.add(String.valueOf(fastq1.getReadLength()));
//        statLine.add(String.valueOf(fastq2.getReadLength()));
//        statLine.add(String.valueOf(consensusCount));
//        statLine.add(String.valueOf(missingCount));
//        statLine.add(String.valueOf(mismatchCount));
//        statLine.add(String.valueOf(modCWithG));
//        statLine.add(String.valueOf(modCWithoutG));
//        statLine.add(String.valueOf(exactMatches));
//        statLine.add(String.valueOf(exactMismatches));
//        statLine.add(String.valueOf(hairpinStart == NO_HAIRPIN ? -1 : hairpinStart + 1));
//        statLine.add(fastq1.getReadString());
//        statLine.add(fastq1.getBaseQualityString());
//        statLine.add(fastq2.getReadString());
//        statLine.add(fastq2.getBaseQualityString());
//        statLine.add(consensusReadForOutput(consensusRead));
//
//        writer.write(statLine.toString());
//        writer.newLine();
//    }

//    private static String consensusReadForOutput(final String consensusRead)
//    {
//        StringBuilder output = new StringBuilder();
//        for(int i = 0; i < consensusRead.length(); i++)
//        {
//            char c = consensusRead.charAt(i);
//            if(c == 'X')
//            {
//                output.append('N');
//                continue;
//            }
//
//            if(c == MODC_BASE)
//            {
//                output.append('X');
//                continue;
//            }
//
//            output.append(c);
//        }
//
//        return output.toString();
//    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        FastqBimodalCollapseConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        FastqBimodalCollapseConfig config = new FastqBimodalCollapseConfig(configBuilder);
        FastqBimodalCollapse fastqBimodalCollapse = new FastqBimodalCollapse(config);
        fastqBimodalCollapse.run();
    }
}
