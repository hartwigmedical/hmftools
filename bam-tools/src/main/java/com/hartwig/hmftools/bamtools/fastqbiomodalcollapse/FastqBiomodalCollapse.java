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

    private static List<BaseQualPair> fastqToSeq(final FastqRecord fastq)
    {
        List<BaseQualPair> seq = Lists.newArrayList();
        for(int i = 0; i < fastq.getReadLength(); i++)
        {
            seq.add(new BaseQualPair(fastq.getReadString().charAt(i), fastq.getBaseQualities()[i]));
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
                return new BaseQualPair(MISSING_BASE, 0);
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
                return new BaseQualPair(MISSING_BASE, 0);
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

    private static String getConsensusRead(final FastqRecord trimmedFastq1, final FastqRecord trimmedFastq2)
    {
        List<BaseQualPair> seq1 = fastqToSeq(trimmedFastq1);
        List<BaseQualPair> seq2 = fastqToSeq(trimmedFastq2);

        StringBuilder consensusRead = new StringBuilder();
        for(int i = 0; i < min(seq1.size(), seq2.size()); i++)
        {
            consensusRead.append(consensusBaseQualPair(seq1.get(i), seq2.get(i)).Base);
        }

        return consensusRead.toString();
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

    private void processFastqPair(final BufferedWriter writer, final FastqRecord fastq1, final FastqRecord fastq2) throws IOException
    {
        if(mFastqPairProcessedCount % 100 == 0)
        {
            System.out.println(mFastqPairProcessedCount);
        }

        mFastqPairProcessedCount++;

        List<BaseQualPair> seq1 = fastqToSeq(fastq1);
        List<BaseQualPair> seqRevComp2 = fastqToSeq(revCompFastq(fastq2));

        List<Pair<BaseQualPair, BaseQualPair>> fullFragmentSeq = ALIGNER.align(seq1, seqRevComp2, FastqBiomodalCollapse::getExactScore, OPEN_GAP_PENALTY, EXTEND_GAP_PENALTY, true, false, false, true);

        List<BaseQualPair> fullFragmentConsensusSeq = collapseAlignedSeqExact(fullFragmentSeq);

        String fullFragmentCigar = getCigar(fullFragmentSeq);

        HairpinInfo fullFragmentHairpin = findHairpin(seqReadString(fullFragmentConsensusSeq), FORWARD_HAIRPIN);
        if(fullFragmentHairpin != null && fullFragmentHairpin.PrefixMatchLength != -1)
        {
            fullFragmentHairpin = null;
        }

        List<BaseQualPair> read1 = Lists.newArrayList();
        List<BaseQualPair> read2 = Lists.newArrayList();
        List<BaseQualPair> read2_b = Lists.newArrayList();
        List<BaseQualPair> read1_b = Lists.newArrayList();
        List<Pair<BaseQualPair, BaseQualPair>> alignedSeq = null;
        List<BaseQualPair> consensusSeq = null;
        String cigar = null;
        List<Pair<BaseQualPair, BaseQualPair>> alignedSeqRevComp = null;
        List<BaseQualPair> consensusSeqRevComp = null;
        String cigarRevComp = null;
        List<Pair<BaseQualPair, BaseQualPair>> finalAlignedSeq = null;
        List<BaseQualPair> finalConsensusSeq = null;
        String finalCigar = null;
        if(fullFragmentHairpin != null)
        {
            int hairpinStart = fullFragmentHairpin.StartIndex;
            int hairpinEnd = fullFragmentHairpin.StartIndex + FORWARD_HAIRPIN.length() - 1;
            for(int i = 0; i < fullFragmentSeq.size(); i++)
            {
                if(i < hairpinStart)
                {
                    if(fullFragmentSeq.get(i).getLeft() != null)
                        read1.add(fullFragmentSeq.get(i).getLeft());

                    if(fullFragmentSeq.get(i).getRight() != null)
                        read2_b.add(fullFragmentSeq.get(i).getRight());

                    continue;
                }

                if(i > hairpinEnd)
                {
                    if(fullFragmentSeq.get(i).getLeft() != null)
                        read1_b.add(fullFragmentSeq.get(i).getLeft().complementBase());

                    if(fullFragmentSeq.get(i).getRight() != null)
                        read2.add(fullFragmentSeq.get(i).getRight().complementBase());
                }
            }

            if(read1.isEmpty() || read2.isEmpty())
            {
                fullFragmentHairpin = null;
            }
            else
            {
                Collections.reverse(read1_b);
                Collections.reverse(read2);

                alignedSeq = ALIGNER.align(read1, read2, FastqBiomodalCollapse::getModCScore, OPEN_GAP_PENALTY, EXTEND_GAP_PENALTY, false, false, false, false);
                consensusSeq = collapseAlignedSeq(alignedSeq);
                cigar = getCigar(alignedSeq);

                alignedSeqRevComp = ALIGNER.align(read2_b, read1_b, FastqBiomodalCollapse::getModCScore, OPEN_GAP_PENALTY, EXTEND_GAP_PENALTY, true, true, false, false);
                consensusSeqRevComp = collapseAlignedSeq(alignedSeqRevComp);
                cigarRevComp = getCigar(alignedSeqRevComp);

                finalAlignedSeq = ALIGNER.align(consensusSeq, consensusSeqRevComp, FastqBiomodalCollapse::getExactScore, OPEN_GAP_PENALTY, EXTEND_GAP_PENALTY, true, false, false, true);
                finalConsensusSeq = collapseAlignedSeqExact(finalAlignedSeq);
                finalCigar = getCigar(finalAlignedSeq);
            }
        }

        if(fullFragmentHairpin == null)
        {
            read1 = fastqToSeq(revCompFastq(fastq1));
            read2 = fastqToSeq(revCompFastq(fastq2));
            read2_b = fastqToSeq(revCompFastq(fastq2));
            read1_b = fastqToSeq(revCompFastq(fastq1));

            alignedSeq = ALIGNER.align(read1, read2, FastqBiomodalCollapse::getModCScore, OPEN_GAP_PENALTY, EXTEND_GAP_PENALTY, false, false, true, true);
            consensusSeq = collapseAlignedSeq(alignedSeq);
            cigar = getCigar(alignedSeq);

            alignedSeqRevComp = ALIGNER.align(read2_b, read1_b, FastqBiomodalCollapse::getModCScore, OPEN_GAP_PENALTY, EXTEND_GAP_PENALTY, true, true, false, false);
            consensusSeqRevComp = collapseAlignedSeq(alignedSeqRevComp);
            cigarRevComp = getCigar(alignedSeqRevComp);

            finalAlignedSeq = ALIGNER.align(consensusSeq, consensusSeqRevComp, FastqBiomodalCollapse::getExactScore, OPEN_GAP_PENALTY, EXTEND_GAP_PENALTY, true, false, false, true);
            finalConsensusSeq = collapseAlignedSeqExact(finalAlignedSeq);
            finalCigar = getCigar(finalAlignedSeq);
        }

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

        statLine.add(firstOfAlignmentReadString(fullFragmentSeq));
        statLine.add(firstOfAlignmentBaseQualString(fullFragmentSeq));
        statLine.add(secondOfAlignmentReadString(fullFragmentSeq));
        statLine.add(secondOfAlignmentBaseQualString(fullFragmentSeq));
        statLine.add(fullFragmentCigar);
        statLine.add(consensusReadForOutput(seqReadString(fullFragmentConsensusSeq)));
        statLine.add(seqBaseQualString(fullFragmentConsensusSeq));

        fullFragmentHairpin = fullFragmentHairpin != null ? fullFragmentHairpin : new HairpinInfo(-2, -1, -1);
        statLine.add(String.valueOf(fullFragmentHairpin.StartIndex + 1));
        statLine.add(String.valueOf(fullFragmentHairpin.MatchCount));

        statLine.add(firstOfAlignmentReadString(alignedSeq));
        statLine.add(firstOfAlignmentBaseQualString(alignedSeq));
        statLine.add(secondOfAlignmentReadString(alignedSeq));
        statLine.add(secondOfAlignmentBaseQualString(alignedSeq));
        statLine.add(firstOfAlignmentReadString(alignedSeqRevComp));
        statLine.add(firstOfAlignmentBaseQualString(alignedSeqRevComp));
        statLine.add(secondOfAlignmentReadString(alignedSeqRevComp));
        statLine.add(secondOfAlignmentBaseQualString(alignedSeqRevComp));
        statLine.add(cigar);
        statLine.add(consensusReadForOutput(firstOfAlignmentReadString(finalAlignedSeq)));
        statLine.add(firstOfAlignmentBaseQualString(finalAlignedSeq));
        statLine.add(cigarRevComp);
        statLine.add(consensusReadForOutput(secondOfAlignmentReadString(finalAlignedSeq)));
        statLine.add(secondOfAlignmentBaseQualString(finalAlignedSeq));
        statLine.add(finalCigar);
        statLine.add(consensusReadForOutput(seqReadString(finalConsensusSeq)));
        statLine.add(seqBaseQualString(finalConsensusSeq));

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
            "full_aligned_read1",
            "full_aligned_qual1",
            "full_aligned_read2_rev_comp",
            "full_aligned_qual2_rev_comp",
            "full_fragment_cigar",
            "full_fragment_consensus_read",
            "full_fragment_consensus_qual",
            "full_fragment_consensus_hairpin_start_pos",
            "full_fragment_consensus_hairpin_8mer_match_count",
            "read1",
            "qual1",
            "read2",
            "qual2",
            "read2_b",
            "qual2_b",
            "read1_b",
            "qual1_b",
            "cigar",
            "consensus_read",
            "consensus_qual",
            "cigar_rc",
            "consensus_rc_read",
            "consensus_rc_qual",
            "final_cigar",
            "final_consensus_read",
            "final_consensus_qual"
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
            char base1 = fastq1.getReadString().charAt(i);
            int qual1 = fastq1.getBaseQualities()[i];
            char base2 = fastq2.getReadString().charAt(i);
            int qual2 = fastq2.getBaseQualities()[i];
            if(base1 == INS_BASE || base2 == INS_BASE)
            {
                indelCount++;
                continue;
            }
            else if(base1 == MISSING_BASE || base2 == MISSING_BASE)
            {
                missingCount++;
                continue;
            }

            char consensusBase = base2 == 'G' ? MISMATCH_BASE : getConsensusBase(base1, base2);

            if(baseIndex(consensusBase) != -1)
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
