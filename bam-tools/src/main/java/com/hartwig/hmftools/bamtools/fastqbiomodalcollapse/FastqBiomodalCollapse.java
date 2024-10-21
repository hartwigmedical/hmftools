package com.hartwig.hmftools.bamtools.fastqbiomodalcollapse;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASES;
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
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.amber.BaseDepthData;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.commons.compress.utils.Lists;
import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.fastq.FastqRecord;

public class FastqBiomodalCollapse
{
    public static class BaseQualPair
    {
        public final char Base;
        // TODO
        public final int Qual;

        public BaseQualPair(char base, int qual)
        {
            Base = base;
            Qual = qual;
        }

        public BaseQualPair complementBase()
        {
            return new BaseQualPair(swapDnaBase(Base), Qual);
        }

        public char qualChar()
        {
            return (char) (Qual + 33);
        }
    }

    private static final char MISSING_BASE = 'N';
    private static final char MISMATCH_BASE = 'X';
    private static final char MODC_BASE = 'c';
    private static final char INS_BASE = '_';
    private static final char ZERO_QUAL = (char) 33;

    private static final String FORWARD_HAIRPIN = "AATGACGATGCGTTCGAGCATCGTTATT";
    private static final String REVERSE_HAIRPIN = "AATAACGATGCTCGAACGCATCGTCATT";

    private static final int OPEN_GAP_PENALTY = 6;
    private static final int EXTEND_GAP_PENALTY = 1;
    private static final NeedlemanWunschAligner<BaseQualPair> ALIGNER = new NeedlemanWunschAligner<>();

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

        // TODO:
        Random random = new Random(0);

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

                // TODO:
//                if(random.nextDouble() <= 0.1)
//                {
//                    processFastqPair(statWriter, fastq1, fastq2);
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
        BT_LOGGER.info("FastqBimodalCollapse complete, mins({}) fastqPairReadCount({}) fastqPairProcessedCount({})", runTimeMinsStr(startTimeMs), mFastqPairReadCount, mFastqPairProcessedCount);
    }

    private static char getConsensusBase(char base1, char base2)
    {
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

    private static int getModCScore(final BaseQualPair base1, final BaseQualPair base2)
    {
        if(base1.Base == MISSING_BASE || base2.Base == MISSING_BASE)
        {
            return 0;
        }

        if(getConsensusBase(base1.Base, base2.Base) != MISMATCH_BASE)
        {
            return 1;
        }

        return -1;
    }

    private static int getModCSwapCompScore(final BaseQualPair base1, final BaseQualPair base2)
    {
        return getModCScore(base2.complementBase(), base1.complementBase());
    }

    public static int getExactScore(final BaseQualPair base1, final BaseQualPair base2)
    {
        if(base1.Base == MISSING_BASE || base2.Base == MISSING_BASE)
        {
            return 0;
        }

        if(base1.Base == base2.Base)
        {
            return 1;
        }

        return -1;
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

        int maxShift = max(read1.length(), read2.length());
        Float bestMismatchProp = null;
        int bestRead1Shift = 0;
        int bestMismatchCount = Integer.MAX_VALUE;
        for(int read1Shift = -maxShift; read1Shift <= min(maxShift, -50); read1Shift++)
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

            if(totalCount < 5)
            {
                continue;
            }

            if(mismatchCount >= 0.2 * totalCount)
            {
                continue;
            }

            float mismatchProp = 1.0f*mismatchCount/totalCount;

            if(bestMismatchProp == null || mismatchProp < bestMismatchProp)
            {
                bestMismatchProp = mismatchProp;
                bestRead1Shift = read1Shift;
                bestMismatchCount = mismatchCount;
            }
        }

        if(bestMismatchProp == null)
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

    private static List<BaseQualPair> fastqToSeq(final FastqRecord fastq, int start, int end)
    {
        List<BaseQualPair> seq = Lists.newArrayList();
        for(int i = start; i <= end; i++)
        {
            seq.add(new BaseQualPair(fastq.getReadString().charAt(i), fastq.getBaseQualities()[i]));
        }

        return seq;
    }

    private static List<BaseQualPair> fastqToSeq(final FastqRecord fastq)
    {
        return fastqToSeq(fastq, 0, fastq.getReadLength() - 1);
    }

    private static List<BaseQualPair> stringToSeq(final String s)
    {
        List<BaseQualPair> seq = Lists.newArrayList();
        for(int i = 0; i < s.length(); i++)
        {
            seq.add(new BaseQualPair(s.charAt(i), 0));
        }

        return seq;
    }

    private static FastqRecord revCompFastq(final FastqRecord fastq)
    {
        StringBuilder bases = new StringBuilder(fastq.getReadString());
        bases = bases.reverse();
        StringBuilder quals = new StringBuilder(fastq.getBaseQualityString());
        quals = quals.reverse();

        for(int i = 0; i < bases.length(); i++)
        {
            char base = bases.charAt(i);
            bases.setCharAt(i, swapDnaBase(base));
        }

        return new FastqRecord(fastq.getReadName(), bases.toString(), fastq.getBaseQualityHeader(), quals.toString());
    }

    private static List<BaseQualPair> revCompSeq(final List<BaseQualPair> seq)
    {
        List<BaseQualPair> output = Lists.newArrayList();
        for(int i = seq.size() - 1; i >= 0; i--)
        {
            BaseQualPair base = seq.get(i);
            output.add(base.complementBase());
        }

        return output;
    }

    private static BaseQualPair consensusBaseQualPair(final BaseQualPair base1, final BaseQualPair base2)
    {
        boolean base1Missing = base1 == null || base1.Base == MISSING_BASE;
        boolean base2Missing = base2 == null || base2.Base == MISSING_BASE;

        if(base1Missing && base2Missing)
        {
            return new BaseQualPair(MISSING_BASE, 0);
        }

        if(base1Missing)
        {
            if(base2.Base == 'G')
            {
                return base2;
            }

            Set<Character> consensusBases = Sets.newHashSet();
            for(char b1 : DNA_BASES)
            {
                char consensusBase = getConsensusBase(b1, base2.Base);
                if(consensusBase != MISMATCH_BASE)
                {
                    consensusBases.add(consensusBase);
                }
            }

            if(consensusBases.contains('C'))
            {
                consensusBases.remove(MODC_BASE);
            }

            if(consensusBases.size() == 1)
            {
                return new BaseQualPair(consensusBases.stream().findFirst().orElse(null), base2.Base == 'G' ? base2.Qual / 2 : base2.Qual);
            }

            return new BaseQualPair(MISSING_BASE, 0);
        }

        if(base2Missing)
        {
            Set<Character> consensusBases = Sets.newHashSet();
            for(char b2 : DNA_BASES)
            {
                if(b2 == 'G')
                {
                    continue;
                }

                char consensusBase = getConsensusBase(base1.Base, b2);
                if(consensusBase != MISMATCH_BASE)
                {
                    consensusBases.add(consensusBase);
                }
            }

            if(consensusBases.contains('C'))
            {
                consensusBases.remove(MODC_BASE);
            }

            if(consensusBases.size() == 1)
            {
                return new BaseQualPair(consensusBases.stream().findFirst().orElse(null), base1.Qual);
            }

            return new BaseQualPair(MISSING_BASE, 0);
        }

        char consensusBase = getConsensusBase(base1.Base, base2.Base);
        if(consensusBase != MISMATCH_BASE)
        {
            if(base2.Base == 'G')
            {
                return new BaseQualPair('G', base1.Qual/2);
            }

            return new BaseQualPair(consensusBase, base2.Base == 'G' ? max(base1.Qual, base2.Qual) / 2 : max(base1.Qual, base2.Qual));
        }

        if(base1.Qual > base2.Qual)
        {
            Set<Character> consensusBases = Sets.newHashSet();
            for(char b2 : DNA_BASES)
            {
                if(b2 == 'G')
                {
                    continue;
                }

                char consensusBase_ = getConsensusBase(base1.Base, b2);
                if(consensusBase_ != MISMATCH_BASE)
                {
                    consensusBases.add(consensusBase_);
                }
            }

            if(consensusBases.contains('C'))
            {
                consensusBases.remove(MODC_BASE);
            }

            if(consensusBases.size() == 1)
            {
                return new BaseQualPair(consensusBases.stream().findFirst().orElse(null), base1.Qual - base2.Qual);
            }

            return new BaseQualPair(MISSING_BASE, 0);
        }

        if(base2.Qual > base1.Qual)
        {
            if(base2.Base == 'G')
            {
                return base2;
            }

            Set<Character> consensusBases = Sets.newHashSet();
            for(char b1 : DNA_BASES)
            {
                char consensusBase_ = getConsensusBase(b1, base2.Base);
                if(consensusBase_ != MISMATCH_BASE)
                {
                    consensusBases.add(consensusBase_);
                }
            }

            if(consensusBases.contains('C'))
            {
                consensusBases.remove(MODC_BASE);
            }

            if(consensusBases.size() == 1)
            {
                return new BaseQualPair(consensusBases.stream().findFirst().orElse(null), base2.Base == 'G' ? (base2.Qual - base1.Qual)/2 : base2.Qual - base1.Qual);
            }

            return new BaseQualPair(MISSING_BASE, 0);
        }

        return new BaseQualPair(MISSING_BASE, 0);
    }

    private static BaseQualPair consensusBaseQualPairExact(final BaseQualPair base1, final BaseQualPair base2)
    {
        boolean base1Missing = base1 == null || base1.Base == MISSING_BASE;
        boolean base2Missing = base2 == null || base2.Base == MISSING_BASE;

        if(base1Missing && base2Missing)
        {
            return new BaseQualPair(MISSING_BASE, 0);
        }

        if(base1Missing)
        {
            return base2;
        }

        if(base2Missing)
        {
            return base1;
        }

        if(base1.Base == base2.Base)
        {
            return new BaseQualPair(base1.Base, max(base1.Qual, base2.Qual));
        }

        if(base1.Qual > base2.Qual)
        {
            return new BaseQualPair(base1.Base, base1.Qual - base2.Qual);
        }

        if(base2.Qual > base1.Qual)
        {
            return new BaseQualPair(base2.Base, base2.Qual - base1.Qual);
        }

        return new BaseQualPair(MISSING_BASE, 0);
    }

    private static List<BaseQualPair> collapseAlignedSeq(final List<Pair<BaseQualPair, BaseQualPair>> alignedSeq)
    {
        return alignedSeq.stream().map(x -> consensusBaseQualPair(x.getLeft(), x.getRight())).collect(Collectors.toList());
    }

    private static List<BaseQualPair> collapseAlignedSeqExact(final List<Pair<BaseQualPair, BaseQualPair>> alignedSeq)
    {
        return alignedSeq.stream().map(x -> consensusBaseQualPairExact(x.getLeft(), x.getRight())).collect(Collectors.toList());
    }

    private static String firstOfAlignmentReadString(final List<Pair<BaseQualPair, BaseQualPair>> alignedSeq)
    {
        return alignedSeq.stream().map(x -> x.getLeft()).map(x -> x == null ? "_" : String.valueOf(x.Base)).collect(Collectors.joining());
    }

    private static String firstOfAlignmentBaseQualString(final List<Pair<BaseQualPair, BaseQualPair>> alignedSeq)
    {
        return alignedSeq.stream().map(x -> x.getLeft()).map(x -> x == null ? String.valueOf(ZERO_QUAL) : String.valueOf(x.qualChar())).collect(Collectors.joining());
    }

    private static String secondOfAlignmentReadString(final List<Pair<BaseQualPair, BaseQualPair>> alignedSeq)
    {
        return alignedSeq.stream().map(x -> x.getRight()).map(x -> x == null ? "_" : String.valueOf(x.Base)).collect(Collectors.joining());
    }

    private static String secondOfAlignmentBaseQualString(final List<Pair<BaseQualPair, BaseQualPair>> alignedSeq)
    {
        return alignedSeq.stream().map(x -> x.getRight()).map(x -> x == null ? String.valueOf(ZERO_QUAL) : String.valueOf(x.qualChar())).collect(Collectors.joining());
    }

    private static String seqReadString(final List<BaseQualPair> seq)
    {
        return seq.stream().map(x -> x == null ? "_" : String.valueOf(x.Base)).collect(Collectors.joining());
    }

    private static String seqBaseQualString(final List<BaseQualPair> seq)
    {
        return seq.stream().map(x -> x == null ? String.valueOf(ZERO_QUAL) : String.valueOf(x.qualChar())).collect(Collectors.joining());
    }

    private static class AlignedHairpinInfo
    {
        public final int FirstHairpinPos;
        public final int LastHairpinPos;
        public final int ExactMatchCount;
        public final int MismatchCount;
        public final int IndelCount;
        public final int SuffixCount;

        public AlignedHairpinInfo(final int firstHairpinPos, final int lastHairpinPos, final int exactMatchCount, final int mismatchCount,
                final int indelCount, final int suffixCount)
        {
            FirstHairpinPos = firstHairpinPos;
            LastHairpinPos = lastHairpinPos;
            ExactMatchCount = exactMatchCount;
            MismatchCount = mismatchCount;
            IndelCount = indelCount;
            SuffixCount = suffixCount;
        }
    }

    private static AlignedHairpinInfo getAlignedHairpinInfo(final List<Pair<BaseQualPair, BaseQualPair>> alignment)
    {
        int matchCount = alignment.stream().mapToInt(x -> x.getLeft() != null && x.getRight() != null ? 1 : 0).sum();
        if(matchCount == 0)
        {
            return null;
        }

        int firstHairpinPos = 0;
        while(alignment.get(firstHairpinPos).getRight() == null)
        {
            firstHairpinPos++;
        }

        int suffixCount = 0;
        while(alignment.get(alignment.size() - 1).getLeft() == null)
        {
            suffixCount++;
            alignment.remove(alignment.size() - 1);
        }

        int lastHairpinPos = alignment.size() - 1;
        while(alignment.get(lastHairpinPos).getRight() == null)
        {
            lastHairpinPos--;
        }

        int exactMatchCount = 0;
        int mismatchCount = 0;
        int indelCount = 0;
        for(int i = firstHairpinPos; i <= lastHairpinPos; i++)
        {
            Pair<BaseQualPair, BaseQualPair> alignedBase = alignment.get(i);
            BaseQualPair base1 = alignedBase.getLeft();
            BaseQualPair base2 = alignedBase.getRight();
            if(base1 == null || base2 == null)
            {
                indelCount++;
                continue;
            }

            if(base1.Base == base2.Base)
            {
                exactMatchCount++;
                continue;
            }

            mismatchCount++;
        }

        return new AlignedHairpinInfo(firstHairpinPos, lastHairpinPos, exactMatchCount, mismatchCount, indelCount, suffixCount);
    }

    private static <T, S> List<Pair<T, S>> zip(List<T> xs, List<S> ys)
    {
        List<Pair<T, S>> output = Lists.newArrayList();
        for(int i = 0; i < max(xs.size(), ys.size()); i++)
        {
            T left = null;
            S right = null;
            if(i < xs.size())
            {
                left = xs.get(i);
            }

            if(i < ys.size())
            {
                right = ys.get(i);
            }

            output.add(Pair.of(left, right));
        }

        return output;
    }

    private void processFastqPair(final BufferedWriter writer, final FastqRecord fastq1_, final FastqRecord fastq2_) throws IOException
    {
        if(mFastqPairProcessedCount % 100 == 0)
        {
            System.out.println(mFastqPairProcessedCount);
        }

        mFastqPairProcessedCount++;

        int trimmedLength = min(fastq1_.getReadLength(), fastq2_.getReadLength());
        FastqRecord fastq1 = new FastqRecord(fastq1_.getReadName(), fastq1_.getReadString().substring(0, trimmedLength), fastq1_.getBaseQualityHeader(), fastq1_.getBaseQualityString().substring(0, trimmedLength));
        FastqRecord fastq2 = new FastqRecord(fastq2_.getReadName(), fastq2_.getReadString().substring(0, trimmedLength), fastq2_.getBaseQualityHeader(), fastq2_.getBaseQualityString().substring(0, trimmedLength));

        List<BaseQualPair> seq1 = fastqToSeq(fastq1);
        List<BaseQualPair> seq2 = fastqToSeq(fastq2);

        List<BaseQualPair> seq1RC = fastqToSeq(revCompFastq(fastq1));
        List<BaseQualPair> seq2RC = fastqToSeq(revCompFastq(fastq2));

        RevCompMatchInfo rcMatch1 = findBestRevCompMatch(fastq1.getReadString(), fastq2.getReadString());
        RevCompMatchInfo rcMatch2 = findBestRevCompMatch(fastq2.getReadString(), fastq1.getReadString());

        HairpinInfo hairpin1 = findHairpin(fastq1.getReadString(), FORWARD_HAIRPIN);
        HairpinInfo hairpin2 = findHairpin(fastq2.getReadString(), REVERSE_HAIRPIN);

        // naive alignment
        List<Pair<BaseQualPair, BaseQualPair>> naiveAlignment = Lists.newArrayList();
        for(int i = 0; i < min(trimmedLength, min(hairpin1 == null ? Integer.MAX_VALUE : hairpin1.StartIndex, hairpin2 == null ? Integer.MAX_VALUE : hairpin2.StartIndex)); i++)
        {
            BaseQualPair base1 = seq1.get(i);
            BaseQualPair base2 = seq2.get(i);
            naiveAlignment.add(Pair.of(base1, base2));
        }

        List<BaseQualPair> naiveConsensus = collapseAlignedSeq(naiveAlignment);
        ReadPairStats naiveStats = getReadPairStats(naiveAlignment);

        // forward naive consensus
        int hairpin1StartIndex = hairpin1 == null ? trimmedLength : hairpin1.StartIndex;
        int hairpin2StartIndex = hairpin2 == null ? trimmedLength : hairpin2.StartIndex;
        int hairpinStartIndex = min(hairpin1StartIndex, hairpin2StartIndex);
        List<BaseQualPair> read1 = seq1.subList(0, hairpinStartIndex);
        List<BaseQualPair> read2 = seq2.subList(0, hairpinStartIndex);
        List<BaseQualPair> forwardConsensus = collapseAlignedSeq(zip(read1, read2));

        // reverse naive consensus
        List<BaseQualPair> reverseConsensus = null;
        if(rcMatch1.MismatchCount != -1)
        {
            int rcStartIndex = -rcMatch1.Read1Shift;
            int rcLength = min(trimmedLength, hairpinStartIndex - rcStartIndex);
            if(rcLength > 0)
            {
                List<BaseQualPair> read1RC = seq1RC.subList(0, rcLength);
                List<BaseQualPair> read2RC = seq2RC.subList(0, rcLength);
                reverseConsensus = collapseAlignedSeq(zip(read2RC, read1RC));
            }
        }

        // forward aligned consensus
        List<Pair<BaseQualPair, BaseQualPair>> forwardAlignment = null;
        read1 = hairpin1 != null ? seq1.subList(0, hairpin1.StartIndex) : seq1;
        read2 = hairpin2 != null ? seq2.subList(0, hairpin2.StartIndex) : seq2;
        if(read1.size() == read2.size())
        {
            forwardAlignment = ALIGNER.align(read1, read2, FastqBiomodalCollapse::getModCScore, OPEN_GAP_PENALTY, EXTEND_GAP_PENALTY, false, false, false, false);
        }
        else if(read1.size() < read2.size())
        {
            forwardAlignment = ALIGNER.align(read1, read2, FastqBiomodalCollapse::getModCScore, OPEN_GAP_PENALTY, EXTEND_GAP_PENALTY, false, false, false, true);
            while(!forwardAlignment.isEmpty() && forwardAlignment.get(forwardAlignment.size() - 1).getLeft() == null)
            {
                forwardAlignment.remove(forwardAlignment.size() - 1);
            }
        }
        else
        {
            forwardAlignment = ALIGNER.align(read1, read2, FastqBiomodalCollapse::getModCScore, OPEN_GAP_PENALTY, EXTEND_GAP_PENALTY, false, false, true, false);
            while(!forwardAlignment.isEmpty() && forwardAlignment.get(forwardAlignment.size() - 1).getRight() == null)
            {
                forwardAlignment.remove(forwardAlignment.size() - 1);
            }
        }

        List<BaseQualPair> forwardAlignmentConsensus = collapseAlignedSeq(forwardAlignment);
        ReadPairStats indelStats = getReadPairStats(forwardAlignment);

        int forwardMatchCount = 0;
        int forwardInsert1Count = 0;
        int forwardInsert2Count = 0;
        for(int i = 0; i < forwardAlignment.size(); i++)
        {
            BaseQualPair base1 = forwardAlignment.get(i).getLeft();
            BaseQualPair base2 = forwardAlignment.get(i).getRight();
            if(base1 != null && base2 != null)
            {
                forwardMatchCount++;
            }
            else if(base1 != null)
            {
                forwardInsert1Count++;
            }
            else
            {
                forwardInsert2Count++;
            }
        }

        // reverse aligned consensus

        // write output
        StringJoiner statLine = new StringJoiner(STAT_DELIMITER);
        statLine.add(fastq1_.getReadName());
        statLine.add(String.valueOf(fastq1_.getReadLength()));
        statLine.add(String.valueOf(fastq2_.getReadLength()));
        statLine.add(fastq1.getReadString());
        statLine.add(fastq1.getBaseQualityString());
        statLine.add(fastq2.getReadString());
        statLine.add(fastq2.getBaseQualityString());

        hairpin1 = hairpin1 != null ? hairpin1 : new HairpinInfo(-2, -1, -1);
        statLine.add(String.valueOf(hairpin1.StartIndex + 1));
        statLine.add(String.valueOf(hairpin1.MatchCount));
        statLine.add(String.valueOf(hairpin1.PrefixMatchLength));

        hairpin2 = hairpin2 != null ? hairpin2 : new HairpinInfo(-2, -1, -1);
        statLine.add(String.valueOf(hairpin2.StartIndex + 1));
        statLine.add(String.valueOf(hairpin2.MatchCount));
        statLine.add(String.valueOf(hairpin2.PrefixMatchLength));

        statLine.add(String.valueOf(rcMatch1.Read1Shift));
        statLine.add(String.valueOf(rcMatch1.MismatchCount));
        statLine.add(String.valueOf(rcMatch2.Read1Shift));
        statLine.add(String.valueOf(rcMatch2.MismatchCount));

        statLine.add(consensusReadForOutput(seqReadString(naiveConsensus)));
        statLine.add(seqBaseQualString(naiveConsensus));

        statLine.add(String.valueOf(naiveStats.MissingCount));
        statLine.add(String.valueOf(naiveStats.MismatchCount));
        statLine.add(String.valueOf(naiveStats.HighQualMismatchCountGG));
        statLine.add(String.valueOf(naiveStats.HighQualMismatchCountOther));
        statLine.add(String.valueOf(naiveStats.LowQualUnambiguousCount));
        statLine.add(String.valueOf(naiveStats.LowQualAmbiguousCount));
        statLine.add(String.valueOf(naiveStats.ModCGCount));
        statLine.add(String.valueOf(naiveStats.ModCOtherCount));

        statLine.add(revCompFastq(fastq1).getReadString());
        statLine.add(revCompFastq(fastq1).getBaseQualityString());
        statLine.add(revCompFastq(fastq2).getReadString());
        statLine.add(revCompFastq(fastq2).getBaseQualityString());

        statLine.add(forwardConsensus == null ? "-" : consensusReadForOutput(seqReadString(forwardConsensus)));
        statLine.add(forwardConsensus == null ? "-" : seqBaseQualString(forwardConsensus));

        statLine.add(reverseConsensus == null ? "-" : consensusReadForOutput(seqReadString(reverseConsensus)));
        statLine.add(reverseConsensus == null ? "-" : seqBaseQualString(reverseConsensus));

        statLine.add(consensusReadForOutput(firstOfAlignmentReadString(forwardAlignment)));
        statLine.add(firstOfAlignmentBaseQualString(forwardAlignment));
        statLine.add(consensusReadForOutput(secondOfAlignmentReadString(forwardAlignment)));
        statLine.add(secondOfAlignmentBaseQualString(forwardAlignment));

        statLine.add(getCigar(forwardAlignment));
        statLine.add(String.valueOf(forwardMatchCount));
        statLine.add(String.valueOf(forwardInsert2Count));
        statLine.add(String.valueOf(forwardInsert1Count));
        statLine.add(String.valueOf(forwardInsert1Count + forwardInsert2Count));

        statLine.add(consensusReadForOutput(seqReadString(forwardAlignmentConsensus)));
        statLine.add(seqBaseQualString(forwardAlignmentConsensus));

        statLine.add(String.valueOf(indelStats.MissingCount));
        statLine.add(String.valueOf(indelStats.MismatchCount));
        statLine.add(String.valueOf(indelStats.HighQualMismatchCountGG));
        statLine.add(String.valueOf(indelStats.HighQualMismatchCountOther));
        statLine.add(String.valueOf(indelStats.LowQualUnambiguousCount));
        statLine.add(String.valueOf(indelStats.LowQualAmbiguousCount));
        statLine.add(String.valueOf(indelStats.ModCGCount));
        statLine.add(String.valueOf(indelStats.ModCOtherCount));

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
            "naive_consensus_read",
            "naive_consensus_qual",
            "missing_count",
            "mismatch_count",
            "high_qual_GG_mismatch_count",
            "high_qual_other_mismatch_count",
            "low_qual_unambiguous_count",
            "low_qual_ambiguous_count",
            "methC_G_count",
            "methC_other_count",
            "read1_rc",
            "qual1_rc",
            "read2_rc",
            "qual2_rc",
            "forward_consensus_read",
            "forward_consensus_qual",
            "reverse_consensus_read",
            "reverse_consensus_qual",
            "aligned_read1",
            "aligned_qual1",
            "aligned_read2",
            "aligned_qual2",
            "forward_cigar",
            "forward_match_count",
            "forward_insert1_count",
            "forward_insert2_count",
            "forward_indel_count",
            "forward_consensus_read",
            "forward_consensus_qual",
            "aligned_missing_count",
            "aligned_mismatch_count",
            "aligned_high_qual_GG_mismatch_count",
            "aligned_high_qual_other_mismatch_count",
            "aligned_low_qual_unambiguous_count",
            "aligned_low_qual_ambiguous_count",
            "aligned_methC_G_count",
            "aligned_methC_other_count"
    };

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

    private static ReadPairStats getReadPairStats(List<Pair<BaseQualPair, BaseQualPair>> alignment)
    {
        int missingCount = 0;
        int mismatchCount = 0;
        int highQualMismatchCountGG = 0;
        int highQualMismatchCountOther = 0;
        int lowQualUnambiguousCount = 0;
        int lowQualAmbiguousCount = 0;
        int modCGCount = 0;
        int modCOtherCount = 0;
        int indelCount = 0;
        for(int i = 0; i < alignment.size(); i++)
        {
            BaseQualPair base1 = alignment.get(i).getLeft();
            BaseQualPair base2 = alignment.get(i).getRight();
            if(base1 == null || base2 == null)
            {
                indelCount++;
                continue;
            }
            else if(base1.Base == MISSING_BASE || base2.Base == MISSING_BASE)
            {
                missingCount++;
                continue;
            }

            char consensusBase;
            consensusBase = base2.Base == 'G' ? MISMATCH_BASE : getConsensusBase(base1.Base, base2.Base);

            if(baseIndex(consensusBase) != -1)
            {
            }
            else if(consensusBase == MISMATCH_BASE)
            {
                mismatchCount++;
                if(min(base1.Qual, base2.Qual) > 30)
                {
                    if(base1.Base == 'G' && base2.Base == 'G')
                    {
                        highQualMismatchCountGG++;
                    }
                    else
                    {
                        highQualMismatchCountOther++;
                    }
                }
                else if(base1.Qual > 30 && base1.Base != 'T')
                {
                    lowQualUnambiguousCount++;
                }
                else if(base2.Qual > 30 && base2.Base != 'A')
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
                if(i == alignment.size() - 1)
                {
                    modCOtherCount++;
                    continue;
                }

                BaseQualPair nextBase1 = alignment.get(i + 1).getLeft();
                BaseQualPair nextBase2 = alignment.get(i + 1).getRight();
                if(nextBase1 == null || nextBase2 == null)
                {
                    modCOtherCount++;
                    continue;
                }

                if(consensusBaseQualPair(nextBase1, nextBase2).Base == 'G')
                {
                    modCGCount++;
                    continue;
                }

                modCOtherCount++;
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
