package com.hartwig.hmftools.bamtools.fastqbiomodalcollapse;

import static java.lang.Math.abs;
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
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.commons.compress.utils.Lists;
import org.apache.commons.lang3.NotImplementedException;
import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.fastq.FastqRecord;

public class FastqBiomodalCollapse
{
    public static class BaseQualPair
    {
        public final char Base;
        public final char Qual;

        public BaseQualPair(char base, char qual)
        {
            Base = base;
            Qual = qual;
        }

        public BaseQualPair complementBase()
        {
            return new BaseQualPair(swapDnaBase(Base), Qual);
        }
    }

    private static final char MISSING_BASE = 'N';
    private static final char MISMATCH_BASE = 'X';
    private static final char MODC_BASE = 'c';
    private static final char INS_BASE = '_';
    private static final char ZERO_QUAL = (char) 33;

    private static final int ROLLING_WINDOW_SIZE = 20;
    private static final int MIN_GOOD_ROLLING_COUNT = 18;
    private static final int NO_HAIRPIN = Integer.MAX_VALUE;

    private static final float HAIRPIN_MATCH_TARGET = 0.9f;
    private static final int MIN_HAIRPIN_MATCH_LENGTH = 7;
    private static final String FORWARD_HAIRPIN = "AATGACGATGCGTTCGAGCATCGTTATT";
    private static final String REVERSE_HAIRPIN = "AATAACGATGCTCGAACGCATCGTCATT";

    private static final NeedlemanWunschAligner<BaseQualPair> ALIGNER = new NeedlemanWunschAligner<>(FastqBiomodalCollapse::getScore, FastqBiomodalCollapse::getScoreSwapComp, 6, 1);

    private final FastqBiomodalCollapseConfig mConfig;

    private int mFastqPairReadCount;
    private int mFastqPairProcessedCount;

    public FastqBiomodalCollapse(final FastqBiomodalCollapseConfig config)
    {
        mConfig = config;

        mFastqPairReadCount = 0;
        mFastqPairProcessedCount = 0;
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
            statWriter.write(Arrays.stream(STAT_HEADERS).collect(Collectors.joining(STAT_DELIMITER)));
            statWriter.newLine();

            FastqRecord fastq1 = nextFastqRecord(fastq1Reader);
            FastqRecord fastq2 = nextFastqRecord(fastq2Reader);
            while(fastq1 != null && fastq2 != null)
            {
                mFastqPairReadCount++;

                if(!fastq1.getReadName().equals(fastq2.getReadName()))
                {
                    throw new RuntimeException("Fastq read name mismatch.");
                }

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
        BT_LOGGER.info("FastqBimodalCollapse complete, mins({}) fastqPairReadCount({}) fastqPairProcessedCount({})", runTimeMinsStr(startTimeMs), mFastqPairReadCount, mFastqPairProcessedCount);
    }

    private static char getConsensusBase(char base1, char base2)
    {
        // TODO: Handle N?
        if(baseIndex(base1) == -1)
        {
            throw new RuntimeException("Invalid base1");
        }

        if(baseIndex(base2) == -1)
        {
            throw new RuntimeException("Invalid base2");
        }

        // TODO: This is deanimation error.
        //        if(base2 == 'G')
        //        {
        //            return MISMATCH_BASE;
        //        }

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

    private static int getScore(final BaseQualPair base1, final BaseQualPair base2)
    {
        if(base1.Base == MISSING_BASE || base2.Base == MISSING_BASE || getConsensusBase(base1.Base, base2.Base) != MISMATCH_BASE)
        {
            return 1;
        }

        return -1;
    }

    private static int getScoreSwapComp(final BaseQualPair base1, final BaseQualPair base2)
    {
        return getScore(base2.complementBase(), base1.complementBase());
    }

    private static String getConsensusRead(final FastqRecord fastq1, final FastqRecord fastq2)
    {
        String read1 = fastq1.getReadString();
        String read2 = fastq2.getReadString();
        StringBuilder consensusRead = new StringBuilder();
        int length = min(read1.length(), read2.length());
        for(int i = 0; i < length; i++)
        {
            char base1 = read1.charAt(i);
            char base2 = read2.charAt(i);

            // TODO: Handle missing and indel bases better.
            if(base1 == MISSING_BASE || base1 == INS_BASE)
            {
                consensusRead.append(MISSING_BASE);
                continue;
            }

            if(base2 == MISSING_BASE || base2 == INS_BASE)
            {
                consensusRead.append(MISSING_BASE);
                continue;
            }

            consensusRead.append(getConsensusBase(base1, base2));
        }

        return consensusRead.toString();
    }

    public static class HairpinInfo
    {
        public final int StartIndex;
        public final int MatchCount;
        public final int PrefixMatchLength;

        public HairpinInfo(final int startIndex, int matchCount, int prefixMatchLength)
        {
            StartIndex = startIndex;
            MatchCount = matchCount;
            PrefixMatchLength = prefixMatchLength;
        }
    }

    @Nullable
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
                return new HairpinInfo(bestIndex, bestCount, -1);
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
                return new HairpinInfo(start, -1, hairpinPrefixLength);
            }
        }

        return null;
    }

    @Nullable
    private static HairpinInfo findAlignedHairpin(final List<Pair<BaseQualPair, BaseQualPair>> alignment)
    {
        Map<Integer, Integer> counts = Maps.newHashMap();
        int start = 0;
        int end = start + 7;
        while(end < FORWARD_HAIRPIN.length())
        {
            int index = -1;
            for(int i = 0; i <= alignment.size() - 8; i++)
            {
                boolean allMatch = true;
                for(int j = start; j <= end; j++)
                {
                    Pair<BaseQualPair, BaseQualPair> alignedBase = alignment.get(i + j - start);
                    BaseQualPair base1 = alignedBase.getLeft();
                    BaseQualPair base2 = alignedBase.getRight();

                    if(base1 != null && base1.Base != FORWARD_HAIRPIN.charAt(j))
                    {
                        allMatch = false;
                        break;
                    }

                    if(base2 != null && base2.Base != REVERSE_HAIRPIN.charAt(j))
                    {
                        allMatch = false;
                        break;
                    }
                }

                if(allMatch)
                {
                    index = i;
                    break;
                }
            }

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
                return new HairpinInfo(bestIndex, bestCount, -1);
            }
        }

        // look for matches at the end
        for(int hairpinPrefixLength = 7; hairpinPrefixLength >= 5; hairpinPrefixLength--)
        {
            start = alignment.size() - hairpinPrefixLength;
            if(start < 0)
            {
                continue;
            }

            boolean allMatch = true;
            for(int i = 0; i < hairpinPrefixLength; i++)
            {
                Pair<BaseQualPair, BaseQualPair> alignedBase = alignment.get(start + i);
                BaseQualPair base1 = alignedBase.getLeft();
                BaseQualPair base2 = alignedBase.getRight();

                if(base1 != null && base1.Base != FORWARD_HAIRPIN.charAt(i))
                {
                    allMatch = false;
                    break;
                }

                if(base2 != null && base2.Base != REVERSE_HAIRPIN.charAt(i))
                {
                    allMatch = false;
                    break;
                }
            }

            if(allMatch)
            {
                return new HairpinInfo(start, -1, hairpinPrefixLength);
            }
        }

        return null;
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

            if(mismatchCount >= 0.2 * totalCount)
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

    private static String getCigar(final List<Pair<BaseQualPair, BaseQualPair>> alignedSeq)
    {
        StringBuilder cigarBuilder = new StringBuilder();
        Character currentOp = null;
        int currentLength = 0;
        for(Pair<BaseQualPair, BaseQualPair> alignedBase : alignedSeq)
        {
            char op = 'M';
            if(alignedBase.getLeft() == null)
            {
                op = 'I';
            }
            else if(alignedBase.getRight() == null)
            {
                op = 'D';
            }

            if(currentOp == null)
            {
                currentOp = op;
                currentLength = 1;
                continue;
            }

            if(currentOp == op)
            {
                currentLength++;
                continue;
            }

            cigarBuilder.append(currentLength);
            cigarBuilder.append(currentOp);

            currentOp = op;
            currentLength = 1;
        }

        if(currentLength > 0)
        {
            cigarBuilder.append(currentLength);
            cigarBuilder.append(currentOp);
        }

        return cigarBuilder.toString();
    }

    private List<BaseQualPair> fastqToSeq(final FastqRecord fastq)
    {
        List<BaseQualPair> seq = Lists.newArrayList();
        for(int i = 0; i < fastq.getReadLength(); i++)
        {
            seq.add(new BaseQualPair(fastq.getReadString().charAt(i), fastq.getBaseQualityString().charAt(i)));
        }

        return seq;
    }

    private List<Pair<BaseQualPair, BaseQualPair>> alignFragment(final FastqRecord fastq1, final FastqRecord fastq2)
    {
        List<BaseQualPair> seq1 = fastqToSeq(fastq1);
        List<BaseQualPair> seq2 = fastqToSeq(fastq2);
        return ALIGNER.align(seq1, seq2, true, true);
    }

    private void processFastqPair(final BufferedWriter writer, final FastqRecord fastq1, final FastqRecord fastq2) throws IOException
    {
        if(mFastqPairProcessedCount % 100 == 0)
        {
            System.out.println(mFastqPairProcessedCount);
        }

        mFastqPairProcessedCount++;

        List<Pair<BaseQualPair, BaseQualPair>> alignedSeq = alignFragment(fastq1, fastq2);
        HairpinInfo alignedHairpin = findAlignedHairpin(alignedSeq);

        String alignedRead1 = alignedSeq.stream().map(x -> x.getLeft()).map(x -> x == null ? "_" : String.valueOf(x.Base)).collect(Collectors.joining());
        String alignedQual1 = alignedSeq.stream().map(x -> x.getLeft()).map(x -> x == null ? String.valueOf(ZERO_QUAL) : String.valueOf(x.Qual)).collect(Collectors.joining());
        String alignedRead2 = alignedSeq.stream().map(x -> x.getRight()).map(x -> x == null ? "_" : String.valueOf(x.Base)).collect(Collectors.joining());
        String alignedQual2 = alignedSeq.stream().map(x -> x.getRight()).map(x -> x == null ? String.valueOf(ZERO_QUAL) : String.valueOf(x.Qual)).collect(Collectors.joining());

        String cigar = getCigar(alignedSeq);

        // write output
        StringJoiner statLine = new StringJoiner(STAT_DELIMITER);
        statLine.add(fastq1.getReadName());
        statLine.add(String.valueOf(fastq1.getReadLength()));
        statLine.add(String.valueOf(fastq2.getReadLength()));
        statLine.add(fastq1.getReadString());
        statLine.add(fastq1.getBaseQualityString());
        statLine.add(fastq2.getReadString());
        statLine.add(fastq2.getBaseQualityString());

        HairpinInfo hairpin1 = findHairpin(fastq1.getReadString(), FORWARD_HAIRPIN);
        hairpin1 = hairpin1 != null ? hairpin1 : new HairpinInfo(-2, -1, -1);
        statLine.add(String.valueOf(hairpin1.StartIndex + 1));
        statLine.add(String.valueOf(hairpin1.MatchCount));
        statLine.add(String.valueOf(hairpin1.PrefixMatchLength));

        HairpinInfo hairpin2 = findHairpin(fastq2.getReadString(), REVERSE_HAIRPIN);
        hairpin2 = hairpin2 != null ? hairpin2 : new HairpinInfo(-2, -1, -1);
        statLine.add(String.valueOf(hairpin2.StartIndex + 1));
        statLine.add(String.valueOf(hairpin2.MatchCount));
        statLine.add(String.valueOf(hairpin2.PrefixMatchLength));

        RevCompMatchInfo rcMatch1 = findBestRevCompMatch(fastq1.getReadString(), fastq2.getReadString());
        RevCompMatchInfo rcMatch2 = findBestRevCompMatch(fastq2.getReadString(), fastq1.getReadString());
        statLine.add(String.valueOf(rcMatch1.Read1Shift));
        statLine.add(String.valueOf(rcMatch1.MismatchCount));
        statLine.add(String.valueOf(rcMatch2.Read1Shift));
        statLine.add(String.valueOf(rcMatch2.MismatchCount));

        hairpin1 = findHairpin(fastq1.getReadString(), FORWARD_HAIRPIN);
        hairpin2 = findHairpin(fastq2.getReadString(), REVERSE_HAIRPIN);
        int trimmedEnd1 = hairpin1 == null ? fastq1.getReadLength() - 1 : hairpin1.StartIndex - 1;
        int trimmedEnd2 = hairpin2 == null ? fastq2.getReadLength() - 1 : hairpin2.StartIndex - 1;
        FastqRecord trimmedFastq1 = new FastqRecord(fastq1.getReadName(), fastq1.getReadString().substring(0, trimmedEnd1 + 1), fastq1.getBaseQualityHeader(), fastq1.getBaseQualityString().substring(0, trimmedEnd1 + 1));
        FastqRecord trimmedFastq2 = new FastqRecord(fastq1.getReadName(), fastq2.getReadString().substring(0, trimmedEnd2 + 1), fastq1.getBaseQualityHeader(), fastq2.getBaseQualityString().substring(0, trimmedEnd2 + 1));
        ReadPairStats stats = getReadPairStats(trimmedFastq1, trimmedFastq2);

        statLine.add(consensusReadForOutput(getConsensusRead(trimmedFastq1, trimmedFastq2)));
        statLine.add(String.valueOf(stats.MissingCount));
        statLine.add(String.valueOf(stats.MismatchCount));
        statLine.add(String.valueOf(stats.HighQualMismatchCountGG));
        statLine.add(String.valueOf(stats.HighQualMismatchCountOther));
        statLine.add(String.valueOf(stats.LowQualUnambiguousCount));
        statLine.add(String.valueOf(stats.LowQualAmbiguousCount));
        statLine.add(String.valueOf(stats.ModCGCount));
        statLine.add(String.valueOf(stats.ModCOtherCount));

        statLine.add(alignedRead1);
        statLine.add(alignedQual1);
        statLine.add(alignedRead2);
        statLine.add(alignedQual2);
        statLine.add(cigar);

        alignedHairpin = alignedHairpin != null ? alignedHairpin : new HairpinInfo(-2, -1, -1);
        statLine.add(String.valueOf(alignedHairpin.StartIndex + 1));
        statLine.add(String.valueOf(alignedHairpin.MatchCount));
        statLine.add(String.valueOf(alignedHairpin.PrefixMatchLength));

        writer.write(statLine.toString());
        writer.newLine();
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
            "hairpin1_8mer_match_count",
            "hairpin1_suffix_match_length",
            "hairpin2_start_pos",
            "hairpin2_8mer_match_count",
            "hairpin2_suffix_match_length",
            "rev_comp_read_shift1",
            "rev_comp_mismatch_count1",
            "rev_comp_read_shift2",
            "rev_comp_mismatch_count2",
            "consensus_read",
            "missing_count",
            "mismatch_count",
            "high_qual_GG_mismatch_count",
            "high_qual_other_mismatch_count",
            "low_qual_unambiguous_count",
            "low_qual_ambiguous_count",
            "methC_G_count",
            "methC_other_count",
            "aligned_read1",
            "aligned_qual1",
            "aligned_read2",
            "aligned_qual2",
            "cigar",
            "aligned_hairpin_start_pos",
            "aligned_hairpin_8mer_match_count",
            "aligned_hairpin_suffix_match_length",
//            ,
//            "aligned_consensus_read",
//            "aligned_missing_count",
//            "aligned_mismatch_count",
//            "aligned_high_qual_GG_mismatch_count",
//            "aligned_high_qual_other_mismatch_count",
//            "aligned_low_qual_unambiguous_count",
//            "aligned_low_qual_ambiguous_count",
//            "aligned_methC_G_count",
//            "aligned_methC_other_count",
//            "aligned_indel_count"
    };

//    private void processFastqPair(final BufferedWriter writer, final FastqRecord fastq1, final FastqRecord fastq2) throws IOException
//    {
//
//        List<Pair<BaseQualPair, BaseQualPair>> alignedSeq = alignUpToHairpin(fastq1, fastq2, hairpin1, hairpin2);
//
//        String alignedRead1 =
//                alignedSeq.stream().map(x -> x.getLeft() == null ? "_" : String.valueOf(x.getLeft().Base)).collect(Collectors.joining());
//        String alignedQual1 = alignedSeq.stream()
//                .map(x -> x.getLeft() == null ? String.valueOf(ZERO_QUAL) : String.valueOf(x.getLeft().Qual))
//                .collect(Collectors.joining());
//
//        String alignedRead2 =
//                alignedSeq.stream().map(x -> x.getRight() == null ? "_" : String.valueOf(x.getRight().Base)).collect(Collectors.joining());
//        String alignedQual2 = alignedSeq.stream()
//                .map(x -> x.getRight() == null ? String.valueOf(ZERO_QUAL) : String.valueOf(x.getRight().Qual))
//                .collect(Collectors.joining());
//
//        String cigar = getCigar(alignedSeq);
//
//        FastqRecord alignedFastq1 = new FastqRecord(fastq1.getReadName(), alignedRead1, fastq1.getBaseQualityHeader(), alignedQual1);
//        FastqRecord alignedFastq2 = new FastqRecord(fastq1.getReadName(), alignedRead2, fastq1.getBaseQualityHeader(), alignedQual2);
//
//        // Get stats
//
//
//        ReadPairStats alignedStats = getReadPairStats(alignedFastq1, alignedFastq2);
//
//        statLine.add(alignedRead1);
//        statLine.add(alignedQual1);
//        statLine.add(alignedRead2);
//        statLine.add(alignedQual2);
//        statLine.add(cigar);
//
//        statLine.add(consensusReadForOutput(getConsensusRead(alignedFastq1, alignedFastq2)));
//
//        statLine.add(String.valueOf(alignedStats.MissingCount));
//        statLine.add(String.valueOf(alignedStats.MismatchCount));
//        statLine.add(String.valueOf(alignedStats.HighQualMismatchCountGG));
//        statLine.add(String.valueOf(alignedStats.HighQualMismatchCountOther));
//        statLine.add(String.valueOf(alignedStats.LowQualUnambiguousCount));
//        statLine.add(String.valueOf(alignedStats.LowQualAmbiguousCount));
//        statLine.add(String.valueOf(alignedStats.ModCGCount));
//        statLine.add(String.valueOf(alignedStats.ModCOtherCount));
//        statLine.add(String.valueOf(alignedStats.IndelCount));
//    }

    public static class ReadPairStats
    {
        public final int MissingCount;
        public final int MismatchCount;
        public final int HighQualMismatchCountGG;
        public final int HighQualMismatchCountOther;
        public final int LowQualUnambiguousCount;
        public final int LowQualAmbiguousCount;
        public final int ModCGCount;
        public final int ModCOtherCount;
        public final int IndelCount;

        public ReadPairStats(final int missingCount, final int mismatchCount, final int highQualMismatchCountGG,
                final int highQualMismatchCountOther,
                final int lowQualUnambiguousCount, final int lowQualAmbiguousCount, final int modCGCount, final int modCOtherCount,
                final int indelCount)
        {
            MissingCount = missingCount;
            MismatchCount = mismatchCount;
            HighQualMismatchCountGG = highQualMismatchCountGG;
            HighQualMismatchCountOther = highQualMismatchCountOther;
            LowQualUnambiguousCount = lowQualUnambiguousCount;
            LowQualAmbiguousCount = lowQualAmbiguousCount;
            ModCGCount = modCGCount;
            ModCOtherCount = modCOtherCount;
            IndelCount = indelCount;
        }
    }

    private static ReadPairStats getReadPairStats(final FastqRecord fastq1, final FastqRecord fastq2)
    {
        String consensusRead = getConsensusRead(fastq1, fastq2);

        int missingCount = 0;
        int mismatchCount = 0;
        int highQualMismatchCountGG = 0;
        int highQualMismatchCountOther = 0;
        int lowQualUnambiguousCount = 0;
        int lowQualAmbiguousCount = 0;
        int modCGCount = 0;
        int modCOtherCount = 0;
        int indelCount = 0;
        for(int i = 0; i < consensusRead.length(); i++)
        {
            char consensusBase = consensusRead.charAt(i);
            char base1 = fastq1.getReadString().charAt(i);
            int qual1 = fastq1.getBaseQualities()[i];
            char base2 = fastq2.getReadString().charAt(i);
            int qual2 = fastq2.getBaseQualities()[i];

            if(base1 == INS_BASE || base2 == INS_BASE)
            {
                indelCount++;
            }
            else if(base1 == MISSING_BASE || base2 == MISSING_BASE)
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
                else if(qual1 > 30 && base1 != 'T')
                {
                    lowQualUnambiguousCount++;
                }
                else if(qual2 > 30 && base2 != 'A')
                {
                    lowQualUnambiguousCount++;
                }
                else
                {
                    lowQualAmbiguousCount++;
                }
            }
            else if(consensusBase == MODC_BASE)
            {
                if(i < consensusRead.length() - 1 && consensusRead.charAt(i + 1) == 'G')
                {
                    modCGCount++;
                }
                else
                {
                    modCOtherCount++;
                }
            }
            else
            {
                throw new RuntimeException("Unreachable");
            }
        }

        return new ReadPairStats(missingCount, mismatchCount, highQualMismatchCountGG, highQualMismatchCountOther, lowQualUnambiguousCount, lowQualAmbiguousCount, modCGCount, modCOtherCount, indelCount);
    }

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
        FastqBiomodalCollapseConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        FastqBiomodalCollapseConfig config = new FastqBiomodalCollapseConfig(configBuilder);
        FastqBiomodalCollapse fastqBimodalCollapse = new FastqBiomodalCollapse(config);
        fastqBimodalCollapse.run();
    }
}
