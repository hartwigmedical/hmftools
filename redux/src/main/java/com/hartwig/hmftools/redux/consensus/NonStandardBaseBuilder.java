package com.hartwig.hmftools.redux.consensus;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.bam.CigarUtils.collapseCigarOps;
import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_BYTES;
import static com.hartwig.hmftools.common.codon.Nucleotides.baseIndex;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.DUPLEX_ERROR_QUAL;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.SIMPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;
import static com.hartwig.hmftools.redux.consensus.BaseBuilder.INVALID_POSITION;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.ALIGNMENT_ONLY;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_MATCH;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_MISMATCH;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.EQ;
import static htsjdk.samtools.CigarOperator.H;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;
import static htsjdk.samtools.CigarOperator.X;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.sequencing.SequencingType;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public abstract class NonStandardBaseBuilder
{
    public static byte ANY_BASE = (byte) 'N';

    protected final RefGenome mRefGenome;
    protected int mChromosomeLength;

    public NonStandardBaseBuilder(final RefGenome refGenome)
    {
        mRefGenome = refGenome;
        mChromosomeLength = 0;
    }

    public abstract void buildConsensusRead(final List<SAMRecord> reads, final ConsensusState consensusState, boolean hasIndels);

    public void setChromosomeLength(int chromosomeLength) { mChromosomeLength = chromosomeLength; }

    public static NonStandardBaseBuilder fromSequencingType(final SequencingType sequencingType, final RefGenome refGenome)
    {
        if(sequencingType == SBX)
            return new SbxBuilder(refGenome);

        return null;
    }

    @VisibleForTesting
    public static class ExtendedRefPos implements Comparable<ExtendedRefPos>
    {
        public final int RefPos;
        public final int InsertIndex;

        public ExtendedRefPos(int refPos, int insertIndex)
        {
            RefPos = refPos;
            InsertIndex = insertIndex;
        }

        @Override
        public boolean equals(final Object o)
        {
            if(this == o)
                return true;

            if(!(o instanceof ExtendedRefPos))
                return false;

            final ExtendedRefPos that = (ExtendedRefPos) o;
            return RefPos == that.RefPos && InsertIndex == that.InsertIndex;
        }

        @Override
        public int hashCode()
        {
            return RefPos + 31 * InsertIndex;
        }

        @Override
        public int compareTo(final ExtendedRefPos o)
        {
            if(RefPos != o.RefPos)
                return RefPos - o.RefPos;

            return InsertIndex - o.InsertIndex;
        }
    }

    @VisibleForTesting
    public static class AnnotatedBase
    {
        public final ExtendedRefPos Pos;
        public final byte Base;
        public final byte Qual;
        public final CigarOperator CigarOp;
        public final Set<String> Annotations;

        public AnnotatedBase(final ExtendedRefPos pos, byte base, byte qual, final CigarOperator cigarOp)
        {
            Pos = pos;
            Base = base;
            Qual = qual;
            CigarOp = cigarOp;
            Annotations = Sets.newHashSet();
        }
    }

    private static List<AnnotatedBase> getAnnotatedBases(final SAMRecord read)
    {
        List<AnnotatedBase> annotatedBases = Lists.newArrayList();
        byte[] bases = read.getReadBases();
        byte[] quals = read.getBaseQualities();

        int readIndex = 0;
        int refPos = read.getAlignmentStart() - leftSoftClipLength(read);
        int insertIndex = 0;
        for(CigarElement el : read.getCigar().getCigarElements())
        {
            if(el.getOperator() == H)
                continue;

            boolean isReadBase = el.getOperator().consumesReadBases();
            boolean isRefBase = el.getOperator() == S || el.getOperator().consumesReferenceBases();

            if(isReadBase && isRefBase)
            {
                insertIndex = 0;
                for(int i = 0; i < el.getLength(); i++)
                {
                    ExtendedRefPos pos = new ExtendedRefPos(refPos, insertIndex);
                    annotatedBases.add(new AnnotatedBase(pos, bases[readIndex], quals[readIndex], el.getOperator()));

                    readIndex++;
                    refPos++;
                }

                continue;
            }

            if(isReadBase)
            {
                for(int i = 0; i < el.getLength(); i++)
                {
                    insertIndex++;
                    ExtendedRefPos pos = new ExtendedRefPos(refPos - 1, insertIndex);
                    annotatedBases.add(new AnnotatedBase(pos, bases[readIndex], quals[readIndex], el.getOperator()));
                    readIndex++;
                }

                continue;
            }

            if(isRefBase)
            {
                insertIndex = 0;
                refPos += el.getLength();
                continue;
            }

            throw new RuntimeException("Unreachable");
        }

        return annotatedBases;
    }

    private static Collection<List<AnnotatedBase>> annotatedAndAlignReads(final List<SAMRecord> reads)
    {
        SortedMap<ExtendedRefPos, List<AnnotatedBase>> alignment = Maps.newTreeMap();
        for(SAMRecord read : reads)
        {
            List<AnnotatedBase> annotatedBases = getAnnotatedBases(read);
            for(AnnotatedBase annotatedBase : annotatedBases)
            {
                alignment.computeIfAbsent(annotatedBase.Pos, key -> Lists.newArrayList());
                alignment.get(annotatedBase.Pos).add(annotatedBase);
            }
        }

        return alignment.values();
    }

    private static List<AnnotatedBase> getConsensusBases(
            final Collection<List<AnnotatedBase>> alignment, final Function<List<AnnotatedBase>, AnnotatedBase> determineConsensus)
    {
        return alignment.stream().map(determineConsensus).collect(Collectors.toList());
    }

    private static CigarOperator getConsensusCigarOp(final List<AnnotatedBase> bases)
    {
        if(bases.get(0).Pos.InsertIndex > 0)
            return I;

        for(AnnotatedBase base : bases)
        {
            CigarOperator op = base.CigarOp;
            if(op == M || op == EQ || op == X)
                return M;
        }

        return S;
    }

    private static void updateConsensusState(final List<SAMRecord> reads, final ConsensusState consensusState, boolean hasIndels,
            final List<AnnotatedBase> consensusBases)
    {
        int minUnclippedPosStart = INVALID_POSITION;
        int maxUnclippedPosEnd = 0;
        int minAlignedPosStart = INVALID_POSITION;
        int maxAlignedPosEnd = 0;

        List<Byte> bases = new ArrayList<>(consensusBases.size());
        List<Byte> quals = new ArrayList<>(consensusBases.size());
        List<CigarOperator> cigarOps = new ArrayList<>(consensusBases.size());
        int lastRefPos = INVALID_POSITION;
        int firstNonClippedIndex = INVALID_POSITION;
        int lastNonClippedIndex = INVALID_POSITION;
        int idx = 0;
        for(AnnotatedBase consensusBase : consensusBases)
        {
            if(consensusBase == null)
                continue;

            int refPos = consensusBase.Pos.RefPos;

            if(lastRefPos != INVALID_POSITION && refPos > lastRefPos + 1)
            {
                for(int i = lastRefPos + 1; i < refPos; i++)
                {
                    cigarOps.add(D);
                    idx++;
                }
            }

            bases.add(consensusBase.Base);
            quals.add(consensusBase.Qual);
            cigarOps.add(consensusBase.CigarOp);

            if(firstNonClippedIndex == INVALID_POSITION && !consensusBase.CigarOp.isClipping())
                firstNonClippedIndex = idx;

            if(!consensusBase.CigarOp.isClipping())
                lastNonClippedIndex = max(lastNonClippedIndex, idx);

            if(minUnclippedPosStart == INVALID_POSITION)
                minUnclippedPosStart = refPos;

            if(minAlignedPosStart == INVALID_POSITION && !consensusBase.CigarOp.isClipping())
                minAlignedPosStart = refPos;

            maxUnclippedPosEnd = max(maxUnclippedPosEnd, refPos);
            if(!consensusBase.CigarOp.isClipping())
                maxAlignedPosEnd = max(maxAlignedPosEnd, refPos);

            lastRefPos = refPos;
            idx++;
        }

        for(int i = firstNonClippedIndex; i <= lastNonClippedIndex; i++)
        {
            if(cigarOps.get(i).isClipping())
                cigarOps.set(i, M);
        }

        consensusState.setBoundaries(minUnclippedPosStart, maxUnclippedPosEnd, minAlignedPosStart, maxAlignedPosEnd);
        consensusState.setBaseLength(bases.size());
        for(int i = 0; i < bases.size(); i++)
        {
            consensusState.Bases[i] = bases.get(i);
            consensusState.BaseQualities[i] = quals.get(i);
        }

        consensusState.CigarElements.addAll(collapseCigarOps(cigarOps));

        if(!hasIndels)
        {
            consensusState.setOutcome(ALIGNMENT_ONLY);
        }
        else
        {
            int uniqueCigarCount = reads.stream()
                    .map(SAMRecord::getCigarString)
                    .collect(Collectors.toCollection(Sets::newHashSet)).size();

            consensusState.setOutcome(uniqueCigarCount == 1 ? INDEL_MATCH : INDEL_MISMATCH);
        }
    }

    public static class SbxBuilder extends NonStandardBaseBuilder
    {
        public SbxBuilder(final RefGenome refGenome)
        {
            super(refGenome);
        }

        @Override
        public void buildConsensusRead(final List<SAMRecord> reads, final ConsensusState consensusState, boolean hasIndels)
        {
            String chromosome = reads.get(0).getReferenceName();
            mChromosomeLength = mRefGenome.getChromosomeLength(chromosome);

            Collection<List<AnnotatedBase>> alignment = annotatedAndAlignReads(reads);
            List<AnnotatedBase> consensusBases = getConsensusBases(alignment, records -> determineConsensus(chromosome, records));
            updateConsensusState(reads, consensusState, hasIndels, consensusBases);
        }

        @VisibleForTesting
        @Nullable
        public AnnotatedBase determineConsensus(final String chromosome, final List<AnnotatedBase> bases)
        {
            ExtendedRefPos pos = bases.get(0).Pos;
            int position = pos.RefPos;
            boolean isInsert = pos.InsertIndex > 0;
            if(isInsert || position < 1 || position > mChromosomeLength)
                position = INVALID_POSITION;

            if(bases.size() == 1)
            {
                if(bases.get(0).Qual == (byte) 0 && position != INVALID_POSITION)
                {
                    byte refBase = mRefGenome.getRefBase(chromosome, position);
                    return new AnnotatedBase(pos, refBase, (byte) 1, M);
                }

                return bases.get(0);
            }

            Map<Byte, int[]> baseCountsByQual = Maps.newHashMap();
            for(AnnotatedBase annotatedBase : bases)
            {
                byte base = annotatedBase.Base;
                int baseIdx = baseIndex(base);
                if(baseIdx < 0 || baseIdx >= DNA_BASE_BYTES.length)
                    continue;

                byte qual = annotatedBase.Qual;
                baseCountsByQual.computeIfAbsent(qual, key -> new int[] { 0, 0, 0, 0 });
                int[] baseCounts = baseCountsByQual.get(qual);
                baseCounts[baseIdx]++;
            }

            if(baseCountsByQual.isEmpty())
            {
                if(position != INVALID_POSITION)
                {
                    byte refBase = mRefGenome.getRefBase(chromosome, position);
                    return new AnnotatedBase(pos, refBase, (byte) 1, M);
                }

                if(isInsert)
                    return null;

                return new AnnotatedBase(pos, ANY_BASE, (byte) 0, S);
            }

            if(!baseCountsByQual.containsKey((byte) SIMPLEX_QUAL) && !baseCountsByQual.containsKey((byte) DUPLEX_QUAL))
            {
                if(position != INVALID_POSITION)
                {
                    byte refBase = mRefGenome.getRefBase(chromosome, position);
                    return new AnnotatedBase(pos, refBase, (byte) 1, M);
                }

                CigarOperator cigarOp = isInsert ? I : S;
                return new AnnotatedBase(pos, ANY_BASE, (byte) 0, cigarOp);
            }

            CigarOperator consensusOp = getConsensusCigarOp(bases);
            if(!baseCountsByQual.containsKey((byte) DUPLEX_QUAL))
            {
                int[] simplexCounts = baseCountsByQual.get((byte) SIMPLEX_QUAL);
                int maxIdx = -1;
                int maxCount = -1;
                int totalCount = 0;
                for(int i = 0; i < simplexCounts.length; i++)
                {
                    totalCount += simplexCounts[i];
                    if(simplexCounts[i] > maxCount)
                    {
                        maxCount = simplexCounts[i];
                        maxIdx = i;
                    }
                }

                if(2 * maxCount > totalCount)
                {
                    byte base = DNA_BASE_BYTES[maxIdx];
                    byte qual = (byte) SIMPLEX_QUAL;
                    return new AnnotatedBase(pos, base, qual, consensusOp);
                }

                if(position != INVALID_POSITION)
                {
                    byte refBase = mRefGenome.getRefBase(chromosome, position);
                    return new AnnotatedBase(pos, refBase, (byte) 1, M);
                }

                return new AnnotatedBase(pos, ANY_BASE, (byte) 0, consensusOp);
            }

            int[] duplexCounts = baseCountsByQual.get((byte) DUPLEX_QUAL);
            int maxIdx = -1;
            boolean multipleMax = false;
            int maxCount = -1;
            int totalCount = 0;
            for(int i = 0; i < duplexCounts.length; i++)
            {
                totalCount += duplexCounts[i];
                if(duplexCounts[i] > maxCount)
                {
                    maxCount = duplexCounts[i];
                    maxIdx = i;
                    multipleMax = false;
                }
                else if(duplexCounts[i] == maxCount)
                {
                    multipleMax = true;
                }
            }

            int[] duplexErrorCounts = baseCountsByQual.get((byte) DUPLEX_ERROR_QUAL);
            if(duplexErrorCounts != null)
            {
                for(int duplexErrorCount : duplexErrorCounts)
                    totalCount += duplexErrorCount;
            }

            if(2 * maxCount > totalCount)
            {
                byte base = DNA_BASE_BYTES[maxIdx];
                byte qual = (byte) DUPLEX_QUAL;
                return new AnnotatedBase(pos, base, qual, consensusOp);
            }

            if(multipleMax)
            {
                if(position != INVALID_POSITION)
                {
                    byte refBase = mRefGenome.getRefBase(chromosome, position);
                    return new AnnotatedBase(pos, refBase, (byte) 1, M);
                }

                return new AnnotatedBase(pos, ANY_BASE, (byte) 0, consensusOp);
            }

            byte base = DNA_BASE_BYTES[maxIdx];
            byte qual = (byte) SIMPLEX_QUAL;
            return new AnnotatedBase(pos, base, qual, consensusOp);
        }
    }
}
