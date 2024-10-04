package com.hartwig.hmftools.bamtools.fastqbimodalcollapse;

import static java.lang.Math.min;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.codon.Nucleotides.baseIndex;
import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
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
    private static final int MIN_HAIRPIN_MATCH_LENGTH = 7;
    private static final String FORWARD_HAIRPIN = "AATGACGATGCGTTCGAGCATCGTTATT";
    private static final String REVERSE_HAIRPIN = "AATAACGATGCTCGAACGCATCGTCATT";

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
        public final int MatchCount;
        public final int PrefixMatchLength;

        public HairpinInfo(final int startPos, int matchCount, int prefixMatchLength)
        {
            StartPos = startPos;
            MatchCount = matchCount;
            PrefixMatchLength = prefixMatchLength;
        }
    }

    private static HairpinInfo findHairpin(final String read, final String hairpinSequence)
    {
        Map<Integer, Integer> counts = Maps.newHashMap();
        int start = 0;
        int end = start + 7;
        while(end < hairpinSequence.length())
        {
            String hairpinSubSeq = hairpinSequence.substring(start, end + 1);
            int index = read.indexOf(hairpinSubSeq);
            if(index < 0)
            {
                start++;
                end++;
                continue;
            }

            index -= start;
            int newCount = counts.getOrDefault(index, 0) + 1;
            counts.put(index, newCount);

            start++;
            end++;
        }

        if(!counts.isEmpty())
        {
            int bestIndex = -1;
            int bestCount = -1;
            for(Map.Entry<Integer, Integer> indexCountPair : counts.entrySet())
            {
                int index = indexCountPair.getKey();
                int count = indexCountPair.getValue();
                if(count > bestCount)
                {
                    bestIndex = index;
                    bestCount = count;
                }
            }

            if(bestCount >= 3)
            {
                return new HairpinInfo(bestIndex + 1, bestCount, -1);
            }
        }

        // look for matches at the end
        for(int hairpinPrefixLength = 7; hairpinPrefixLength >= 5; hairpinPrefixLength--)
        {
            start = read.length() - hairpinPrefixLength;
            if(start < 0)
            {
                continue;
            }

            boolean all_matches = true;
            for(int i = 0; i < hairpinPrefixLength; i++)
            {
                char readBase = read.charAt(start + i);
                char hairpinBase = hairpinSequence.charAt(i);
                if(readBase != hairpinBase)
                {
                    all_matches = false;
                    break;
                }
            }

            if(all_matches)
            {
                return new HairpinInfo(start + 1, -1, hairpinPrefixLength);
            }
        }

        return new HairpinInfo(-1, -1, -1);
    }

    public static class RevCompMatchInfo
    {
        public final int Read1Shift;
        public final int MismatchCount;

        public RevCompMatchInfo(int read1Shift, int mismatchCount)
        {
            Read1Shift = read1Shift;
            MismatchCount = mismatchCount;
        }
    }

    private static RevCompMatchInfo findBestRevCompMatch(final String read1, final String read2)
    {
        StringBuilder read2RevCompBuilder = new StringBuilder();
        for(int i = read2.length() - 1; i >= 0; i--)
        {
            read2RevCompBuilder.append(swapDnaBase(read2.charAt(i)));
        }
        String read2RevComp = read2RevCompBuilder.toString();

        int bestMismatchCount = Integer.MAX_VALUE;
        int bestRead1Shift = 0;
        for(int read1Shift = -120; read1Shift <= 120; read1Shift++)
        {
            int i1 = 0;
            int i2 = 0;
            if(read1Shift < 0)
            {
                i1 = -read1Shift;
            }

            if(read1Shift > 0)
            {
                i2 = read1Shift;
            }

            int totalCount = 0;
            int mismatchCount = 0;
            while(i1 < read1.length() && i2 < read2RevComp.length())
            {
                char base1 = read1.charAt(i1);
                char base2 = read2RevComp.charAt(i2);
                if(base1 != base2)
                {
                    mismatchCount++;
                }

                totalCount++;
                i1++;
                i2++;
            }

            if(totalCount == 0)
            {
                continue;
            }

            if(mismatchCount >= 0.2*totalCount)
            {
                continue;
            }

            if(mismatchCount < bestMismatchCount)
            {
                bestMismatchCount = mismatchCount;
                bestRead1Shift = read1Shift;
            }
        }

        if(bestMismatchCount == Integer.MAX_VALUE)
        {
            return new RevCompMatchInfo(0, -1);
        }

        return new RevCompMatchInfo(bestRead1Shift, bestMismatchCount);
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
            "hairpin1_8mer_count",
            "hairpin1_end_match_length",
            "hairpin2_start_pos",
            "hairpin2_8mer_count",
            "hairpin2_end_match_length",
            "rev_comp_read1_shift",
            "rev_comp_mismatch_count",
            "consensus_read",
            "missing_count",
            "mismatch_count",
            "GG_high_qual_mismatch_count",
            "other_high_qual_mismatch_count",
            "methC_count"
    };

    private static void processFastqPair(final BufferedWriter writer, final FastqRecord fastq1, final FastqRecord fastq2) throws IOException
    {
        HairpinInfo hairpin1 = findHairpin(fastq1.getReadString(), FORWARD_HAIRPIN);
        HairpinInfo hairpin2 = findHairpin(fastq2.getReadString(), REVERSE_HAIRPIN);
        RevCompMatchInfo revCompMatchInfo = findBestRevCompMatch(fastq1.getReadString(), fastq2.getReadString());

        String consensusRead = getConsensusRead(fastq1, fastq2);
        int hairpinStart = min(hairpin1.StartPos == -1 ? NO_HAIRPIN : hairpin1.StartPos, hairpin2.StartPos == -1 ? NO_HAIRPIN : hairpin2.StartPos);

        int missingCount = 0;
        int mismatchCount = 0;
        int highQualMismatchCountGG = 0;
        int highQualMismatchCountOther = 0;
        int modCCount = 0;
        for(int i = 0; i < min(consensusRead.length(), hairpinStart); i++)
        {
            char consensusBase = consensusRead.charAt(i);
            char base1 = fastq1.getReadString().charAt(i);
            int qual1 = fastq1.getBaseQualities()[i];
            char base2 = fastq2.getReadString().charAt(i);
            int qual2 = fastq2.getBaseQualities()[i];

            if(base1 == MISSING_BASE || base2 == MISSING_BASE)
            {
                missingCount++;
            }
            else if(baseIndex(consensusBase) != -1)
            {
            }
            else if(consensusBase == MISMATCH_BASE)
            {
                mismatchCount++;
                if(min(qual1, qual2) > 30)
                {
                    if(base1 == 'G' && base2 == 'G')
                    {
                        highQualMismatchCountGG++;
                    }
                    else
                    {
                        highQualMismatchCountOther++;
                    }
                }
            }
            else if(consensusBase == MODC_BASE)
            {
                modCCount++;
            }
            else
            {
                throw new RuntimeException("Unreachable");
            }
        }

        StringJoiner statLine = new StringJoiner(STAT_DELIMITER);
        statLine.add(fastq1.getReadName());
        statLine.add(String.valueOf(fastq1.getReadLength()));
        statLine.add(String.valueOf(fastq2.getReadLength()));
        statLine.add(fastq1.getReadString());
        statLine.add(fastq1.getBaseQualityString());
        statLine.add(fastq2.getReadString());
        statLine.add(fastq2.getBaseQualityString());
        statLine.add(String.valueOf(hairpin1.StartPos));
        statLine.add(String.valueOf(hairpin1.MatchCount));
        statLine.add(String.valueOf(hairpin1.PrefixMatchLength));
        statLine.add(String.valueOf(hairpin2.StartPos));
        statLine.add(String.valueOf(hairpin2.MatchCount));
        statLine.add(String.valueOf(hairpin2.PrefixMatchLength));
        statLine.add(String.valueOf(revCompMatchInfo.Read1Shift));
        statLine.add(String.valueOf(revCompMatchInfo.MismatchCount));
        statLine.add(consensusReadForOutput(consensusRead));
        statLine.add(String.valueOf(missingCount));
        statLine.add(String.valueOf(mismatchCount));
        statLine.add(String.valueOf(highQualMismatchCountGG));
        statLine.add(String.valueOf(highQualMismatchCountOther));
        statLine.add(String.valueOf(modCCount));

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

    private static String consensusReadForOutput(final String consensusRead)
    {
        StringBuilder output = new StringBuilder();
        for(int i = 0; i < consensusRead.length(); i++)
        {
            char c = consensusRead.charAt(i);
            if(c == 'X')
            {
                output.append('N');
                continue;
            }

            if(c == MODC_BASE)
            {
                output.append('X');
                continue;
            }

            output.append(c);
        }

        return output.toString();
    }

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
