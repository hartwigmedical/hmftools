package com.hartwig.hmftools.fastqtools.biomodalcollapse;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.Nucleotides.baseIndex;
import static com.hartwig.hmftools.common.sequencing.BiomodalBamUtils.LOW_QUAL_CUTOFF;
import static com.hartwig.hmftools.fastqtools.biomodalcollapse.BiomodalCollapseUtil.QualCappingOption.CAP_BY_FIRST;
import static com.hartwig.hmftools.fastqtools.biomodalcollapse.BiomodalCollapseUtil.QualCappingOption.CAP_BY_SECOND;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.fastq.FastqRecord;

public class BiomodalCollapseUtil
{
    @Nullable
    public static FastqRecord nextFastqRecord(final BufferedReader reader)
    {
        try
        {
            String rawReadName = reader.readLine();
            if(rawReadName == null)
            {
                return null;
            }

            // ignore tags
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
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }
    }

    public static List<BaseQualPair> fastqToSeq(final FastqRecord fastq, int start, int end)
    {
        List<BaseQualPair> seq = Lists.newArrayList();
        byte[] readBases = fastq.getReadBases();
        byte[] quals = fastq.getBaseQualities();
        for(int i = start; i <= end; i++)
        {
            seq.add(new BaseQualPair(readBases[i], quals[i]));
        }

        return seq;
    }

    public static List<BaseQualPair> fastqToSeq(final FastqRecord fastq)
    {
        return fastqToSeq(fastq, 0, fastq.getReadLength() - 1);
    }

    public static FastqRecord seqToFastq(final String readName, final String qualityHeader, final List<BaseQualPair> clippedSeq)
    {
        byte[] readBases = readBasesFromSeq(clippedSeq);
        byte[] baseQualities = baseQualitiesFromSeq(clippedSeq);
        return new FastqRecord(readName, readBases, qualityHeader, baseQualities);
    }

    public static FastqRecord seqToFastq(final List<BaseQualPair> clippedSeq)
    {
        return seqToFastq("", "", clippedSeq);
    }

    public static byte[] readBasesFromSeq(final List<BaseQualPair> seq)
    {
        byte[] readBases = new byte[seq.size()];
        for(int i = 0; i < seq.size(); i++)
        {
            readBases[i] = seq.get(i).Base;
        }

        return readBases;
    }

    public static byte[] baseQualitiesFromSeq(final List<BaseQualPair> seq)
    {
        byte[] baseQualities = new byte[seq.size()];
        for(int i = 0; i < seq.size(); i++)
        {
            baseQualities[i] = (byte) seq.get(i).Qual;
        }

        return baseQualities;
    }

    public static List<BaseQualPair> reverseComplementSeq(final List<BaseQualPair> seq)
    {
        List<BaseQualPair> output = Lists.newArrayList();
        for(int i = seq.size() - 1; i >= 0; i--)
        {
            output.add(seq.get(i).complementBase());
        }

        return output;
    }

    public static int getCutPoint(int trimmedLength, final HairpinInfo hairpin1, final HairpinInfo hairpin2,
            final ReverseComplementMatchInfo reverseComplementMatchInfo)
    {
        int hairpinCutPoint =
                min(hairpin1 == null ? trimmedLength : hairpin1.StartIndex, hairpin2 == null ? trimmedLength : hairpin2.StartIndex);

        int revCompCutPoint = trimmedLength;
        if(reverseComplementMatchInfo != null)
        {
            if(reverseComplementMatchInfo.Read1Shift <= 0)
            {
                int prefixLength = -reverseComplementMatchInfo.Read1Shift;
                int overlapLength = trimmedLength - prefixLength;
                revCompCutPoint = prefixLength + overlapLength / 2 + (overlapLength % 2 == 0 ? 0 : 1);
            }
            else
            {
                int overlapLength = trimmedLength - reverseComplementMatchInfo.Read1Shift;
                revCompCutPoint = overlapLength / 2 + (overlapLength % 2 == 0 ? 0 : 1);
            }
        }

        return min(hairpinCutPoint, revCompCutPoint);
    }

    public static byte getConsensusBase(byte base1, byte base2)
    {
        if(baseIndex(base1) == -1)
        {
            throw new RuntimeException("Invalid base1");
        }

        if(baseIndex(base2) == -1)
        {
            throw new RuntimeException("Invalid base2");
        }

        // Ignore this mismatch, because it can be caused by a deanimation error
        //        if(base2 == 'G')
        //        {
        //            return MISMATCH_BASE;
        //        }

        if(base1 == base2)
        {
            return base1 == (byte) 'C' ? BiomodalConstants.MODC_BASE : base1;
        }

        if(base1 == (byte) 'G')
        {
            return base2 == (byte) 'A' ? (byte) 'G' : BiomodalConstants.MISMATCH_BASE;
        }

        if(base1 == (byte) 'T')
        {
            return base2 == (byte) 'C' ? (byte) 'C' : BiomodalConstants.MISMATCH_BASE;
        }

        return BiomodalConstants.MISMATCH_BASE;
    }

    public static int getModCScore(final BaseQualPair base1, final BaseQualPair base2)
    {
        if(base1.Base == BiomodalConstants.MISSING_BASE || base2.Base == BiomodalConstants.MISSING_BASE)
        {
            return 0;
        }

        if(getConsensusBase(base1.Base, base2.Base) != BiomodalConstants.MISMATCH_BASE)
        {
            return 1;
        }

        return -1;
    }

    public static int getExactScore(final BaseQualPair base1, final BaseQualPair base2, boolean collapseModC)
    {
        if(base1.Base == BiomodalConstants.MISSING_BASE || base2.Base == BiomodalConstants.MISSING_BASE)
        {
            return 0;
        }

        byte collapsedBase1 = base1.Base == BiomodalConstants.MODC_BASE && collapseModC ? (byte) 'C' : base1.Base;
        byte collapsedBase2 = base2.Base == BiomodalConstants.MODC_BASE && collapseModC ? (byte) 'C' : base2.Base;
        if(collapsedBase1 == collapsedBase2)
        {
            return 1;
        }

        return -1;
    }

    private static Set<Byte> modCConsensusBaseOptionsFromRead1Base(byte base1)
    {
        if(base1 == (byte) 'A')
        {
            return Sets.newHashSet((byte) 'A');
        }

        if(base1 == (byte) 'T')
        {
            return Sets.newHashSet((byte) 'A', (byte) 'C');
        }

        if(base1 == (byte) 'G')
        {
            return Sets.newHashSet((byte) 'G');
        }

        if(base1 == (byte) 'C')
        {
            return Sets.newHashSet(BiomodalConstants.MODC_BASE);
        }

        throw new RuntimeException("Unreachable");
    }

    private static Set<Byte> modCConsensusBaseOptionsFromRead2Base(byte base2)
    {
        if(base2 == (byte) 'A')
        {
            return Sets.newHashSet((byte) 'A', (byte) 'G');
        }

        if(base2 == (byte) 'T')
        {
            return Sets.newHashSet((byte) 'T');
        }

        // Ignore this "error", because it can be caused by a deanimation error
        if(base2 == (byte) 'G')
        {
            return Sets.newHashSet((byte) 'G');
        }

        if(base2 == (byte) 'C')
        {
            // ignore modC
            return Sets.newHashSet((byte) 'C');
        }

        throw new RuntimeException("Unreachable");
    }

    public static BaseQualPair modCConsensusBaseQualPair(final BaseQualPair base1, final BaseQualPair base2)
    {
        boolean base1Missing = base1 == null || base1.Base == BiomodalConstants.MISSING_BASE;
        boolean base2Missing = base2 == null || base2.Base == BiomodalConstants.MISSING_BASE;
        boolean bothBasesExist = !base1Missing && !base2Missing;

        if(base1Missing && base2Missing)
        {
            return new BaseQualPair(BiomodalConstants.MISSING_BASE, 0);
        }

        byte consensusBase = base1Missing || base2Missing ? BiomodalConstants.MISMATCH_BASE : getConsensusBase(base1.Base, base2.Base);
        if(bothBasesExist && consensusBase == BiomodalConstants.MISMATCH_BASE && base1.Qual == base2.Qual)
        {
            return new BaseQualPair(BiomodalConstants.MISSING_BASE, 0);
        }

        if(base1Missing || bothBasesExist && consensusBase == BiomodalConstants.MISMATCH_BASE && base2.Qual > base1.Qual)
        {
            Set<Byte> consensusBases = modCConsensusBaseOptionsFromRead2Base(base2.Base);
            if(consensusBases.size() == 1)
            {
                byte base = consensusBases.stream().findFirst().orElse(null);
                int base2Qual = base2.Base == (byte) 'G' || !bothBasesExist ? base2.Qual / 2 : base2.Qual;
                int qual = max(0, bothBasesExist ? base2Qual - base1.Qual : base2Qual);
                return new BaseQualPair(base, qual);
            }

            return new BaseQualPair(BiomodalConstants.MISSING_BASE, 0);
        }

        if(base2Missing || bothBasesExist && consensusBase == BiomodalConstants.MISMATCH_BASE && base1.Qual > base2.Qual)
        {
            Set<Byte> consensusBases = modCConsensusBaseOptionsFromRead1Base(base1.Base);
            if(consensusBases.size() == 1)
            {
                byte base = consensusBases.stream().findFirst().orElse(null);
                int base1Qual = base1.Base == (byte) 'C' || !bothBasesExist ? base1.Qual / 2 : base1.Qual;
                int qual = max(0, bothBasesExist ? base1Qual - base2.Qual : base1Qual);
                return new BaseQualPair(base, qual);
            }

            return new BaseQualPair(BiomodalConstants.MISSING_BASE, 0);
        }

        if(consensusBase != BiomodalConstants.MISMATCH_BASE)
        {
            if(base2.Base == (byte) 'G')
            {
                return new BaseQualPair((byte) 'G', base1.Qual / 2);
            }

            if(base1.Qual == base2.Qual)
            {
                return new BaseQualPair(consensusBase, base1.Qual);
            }

            Set<Byte> highQualConsensusBases = base1.Qual > base2.Qual
                    ? modCConsensusBaseOptionsFromRead1Base(base1.Base)
                    : modCConsensusBaseOptionsFromRead2Base(base2.Base);
            int qual = highQualConsensusBases.size() > 1 ? min(base1.Qual, base2.Qual) : max(base1.Qual, base2.Qual);
            return new BaseQualPair(consensusBase, qual);
        }

        throw new RuntimeException("Unreachable");
    }

    public enum QualCappingOption
    {
        CAP_BY_FIRST,
        CAP_BY_SECOND,
        NONE
    }

    private static BaseQualPair exactConsensusBaseQualPair(final BaseQualPair base1, final BaseQualPair base2,
            final QualCappingOption qualCapping)
    {
        BaseQualPair consensus = exactConsensusBaseQualPair(base1, base2);
        int qual = consensus.Qual;
        if(qualCapping == CAP_BY_FIRST)
        {
            boolean base1Missing = base1 == null || base1.Base == BiomodalConstants.MISSING_BASE;
            int base1Qual = base1Missing ? 0 : base1.Qual;
            qual = min(qual, base1Qual);
        }
        else if(qualCapping == CAP_BY_SECOND)
        {
            boolean base2Missing = base2 == null || base2.Base == BiomodalConstants.MISSING_BASE;
            int base2Qual = base2Missing ? 0 : base2.Qual;
            qual = min(qual, base2Qual);
        }

        return new BaseQualPair(consensus.Base, qual);
    }

    public static BaseQualPair exactConsensusBaseQualPair(final BaseQualPair base1, final BaseQualPair base2)
    {
        boolean base1Missing = base1 == null || base1.Base == BiomodalConstants.MISSING_BASE;
        boolean base2Missing = base2 == null || base2.Base == BiomodalConstants.MISSING_BASE;

        if(base1Missing && base2Missing)
        {
            return new BaseQualPair(BiomodalConstants.MISSING_BASE, 0);
        }

        if(base1Missing)
        {
            return base2;
        }

        if(base2Missing)
        {
            return base1;
        }

        if(base1.Base == BiomodalConstants.MODC_BASE && base2.Base == BiomodalConstants.MODC_BASE)
        {
            return new BaseQualPair(BiomodalConstants.MODC_BASE, max(base1.Qual, base2.Qual));
        }

        byte collapsedBase1 = base1.Base == BiomodalConstants.MODC_BASE ? (byte) 'C' : base1.Base;
        byte collapsedBase2 = base2.Base == BiomodalConstants.MODC_BASE ? (byte) 'C' : base2.Base;
        if(collapsedBase1 == collapsedBase2)
        {
            return new BaseQualPair(collapsedBase1, max(base1.Qual, base2.Qual));
        }

        if(base1.Qual > base2.Qual)
        {
            return new BaseQualPair(base1.Base, base1.Qual - base2.Qual);
        }

        if(base2.Qual > base1.Qual)
        {
            return new BaseQualPair(base2.Base, base2.Qual - base1.Qual);
        }

        return new BaseQualPair(BiomodalConstants.MISSING_BASE, 0);
    }

    public static List<BaseQualPair> collapseAlignedSeqModC(final List<Pair<BaseQualPair, BaseQualPair>> alignedSeq)
    {
        return alignedSeq.stream().map(x -> modCConsensusBaseQualPair(x.getLeft(), x.getRight())).collect(Collectors.toList());
    }

    public static List<BaseQualPair> collapseAlignedSeqExact(final List<Pair<BaseQualPair, BaseQualPair>> alignedSeq,
            QualCappingOption qualCapping)
    {
        return alignedSeq.stream()
                .map(x -> exactConsensusBaseQualPair(x.getLeft(), x.getRight(), qualCapping))
                .collect(Collectors.toList());
    }

    public static <T, S> List<Pair<T, S>> zip(List<T> xs, List<S> ys)
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

    public static <T, S> List<T> getLeftElements(final List<Pair<T, S>> xs, final T defaultValue)
    {
        return xs.stream().map(x -> x.getLeft()).map(x -> x == null ? defaultValue : x).collect(Collectors.toList());
    }

    public static <T, S> List<S> getRightElements(final List<Pair<T, S>> xs, final S defaultValue)
    {
        return xs.stream().map(x -> x.getRight()).map(x -> x == null ? defaultValue : x).collect(Collectors.toList());
    }

    private static class SuffixClipper
    {
        public final Predicate<BaseQualPair> IsBadBase;
        public final float PropThreshold;

        public SuffixClipper(final Predicate<BaseQualPair> isBadBase, float propThreshold)
        {
            IsBadBase = isBadBase;
            PropThreshold = propThreshold;
        }
    }

    private static int suffixClipLength(final List<BaseQualPair> seq, final List<SuffixClipper> clippers)
    {
        if(seq.isEmpty() || clippers.isEmpty())
        {
            return 0;
        }

        List<List<Float>> cScores = Lists.newArrayList();
        List<List<Float>> minCScores = Lists.newArrayList();
        for(SuffixClipper clipper : clippers)
        {
            float badPenalty = 1.0f / clipper.PropThreshold - 1.0f;
            List<Float> cScore = new ArrayList<>(seq.size() + 1);
            List<Float> minCScore = new ArrayList<>(seq.size());

            cScore.add(0.0f);
            for(int i = 0; i < seq.size(); i++)
            {
                BaseQualPair base = seq.get(seq.size() - 1 - i);
                float baseScore = clipper.IsBadBase.test(base) ? -badPenalty : 1.0f;
                cScore.add(cScore.get(i) + baseScore);
            }

            minCScore.add(cScore.get(seq.size()));
            for(int i = 1; i < seq.size(); i++)
            {
                minCScore.add(min(minCScore.get(i - 1), cScore.get(seq.size() - i)));
            }

            cScores.add(cScore);
            minCScores.add(minCScore);
        }

        for(int i = 0; i < seq.size(); i++)
        {
            boolean allGood = true;
            for(int j = 0; j < cScores.size(); j++)
            {
                if(minCScores.get(j).get(seq.size() - 1 - i) - cScores.get(j).get(i) < 0)
                {
                    allGood = false;
                    break;
                }
            }

            if(allGood)
            {
                return i;
            }
        }

        return seq.size();
    }

    public static Pair<Integer, Integer> qualityTrim(final List<BaseQualPair> seq)
    {
        List<SuffixClipper> suffixClippers = Lists.newArrayList();
        suffixClippers.add(new SuffixClipper(x -> x.Base != BiomodalConstants.MISSING_BASE
                && x.Qual <= LOW_QUAL_CUTOFF, BiomodalConstants.LOW_QUAL_TRIM_PROPORTION_THRESHOLD));
        suffixClippers.add(new SuffixClipper(x -> x.Base == BiomodalConstants.MISSING_BASE, BiomodalConstants.MISSING_BASE_TRIM_PROPORTION_THRESHOLD));

        List<BaseQualPair> trimmedSeq = seq.subList(BiomodalConstants.PREFIX_TRIM_LENGTH, seq.size());

        int suffixTrimLength = suffixClipLength(trimmedSeq, suffixClippers);
        if(suffixTrimLength == trimmedSeq.size())
        {
            return null;
        }

        return Pair.of(BiomodalConstants.PREFIX_TRIM_LENGTH, seq.size() - 1 - suffixTrimLength);
    }

    public static String getCigar(final List<Pair<BaseQualPair, BaseQualPair>> alignedSeq)
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

    public static String consensusReadForStatOutput(final String consensusRead)
    {
        StringBuilder output = new StringBuilder();
        for(int i = 0; i < consensusRead.length(); i++)
        {
            char c = consensusRead.charAt(i);
            if(c == (char) BiomodalConstants.MISMATCH_BASE)
            {
                output.append((char) BiomodalConstants.MISSING_BASE);
                continue;
            }

            if(c == (char) BiomodalConstants.MODC_BASE)
            {
                output.append('X');
                continue;
            }

            output.append(c);
        }

        return output.toString();
    }

    public static String sanatizeQualString(final String qualString)
    {
        StringBuilder output = new StringBuilder(qualString);
        for(int i = 0; i < output.length(); i++)
        {
            // double quotes are problematic for tsv files
            if(output.charAt(i) == '"')
            {
                output.setCharAt(i, (char) ((int) '"' - 1));
            }
        }

        return output.toString();
    }

    private static final int PAIRS_TO_CHECK_FOR_SWAP_DETECTION = 10_000;

    public static boolean isSwappedR1R2(final BufferedReader fastq1Reader, final BufferedReader fastq2Reader)
    {
        int cSeq1Count = 0;
        int gSeq1Count = 0;
        int cSeq2Count = 0;
        int gSeq2Count = 0;

        SynchronizedPairedFastqReader fastqPairReader =
                new SynchronizedPairedFastqReader(fastq1Reader, fastq2Reader, PAIRS_TO_CHECK_FOR_SWAP_DETECTION);

        Pair<FastqRecord, FastqRecord> fastqPair;
        while((fastqPair = fastqPairReader.getNext()) != null)
        {
            FastqRecord fastq1 = fastqPair.getLeft();
            FastqRecord fastq2 = fastqPair.getRight();

            int trimmedLength = min(fastq1.getReadLength(), fastq2.getReadLength());

            List<BaseQualPair> seq1 = fastqToSeq(fastq1, 0, trimmedLength - 1);
            List<BaseQualPair> seq2 = fastqToSeq(fastq2, 0, trimmedLength - 1);
            for(int i = 0; i < trimmedLength; i++)
            {
                byte base1 = seq1.get(i).Base;
                byte base2 = seq2.get(i).Base;
                if(base1 == (byte) 'C')
                    cSeq1Count++;
                else if(base1 == (byte) 'G')
                    gSeq1Count++;

                if(base2 == (byte) 'C')
                    cSeq2Count++;
                else if(base2 == (byte) 'G')
                    gSeq2Count++;
            }
        }

        return cSeq1Count > gSeq1Count && gSeq2Count > cSeq2Count;
    }
}
