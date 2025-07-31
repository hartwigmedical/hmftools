package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_CLIPPING_PENALTY;
import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_GAP_EXTEND_PENALTY;
import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_GAP_OPEN_PENALTY;
import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_MATCH_SCORE;
import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_MISMATCH_PENALTY;
import static com.hartwig.hmftools.common.bam.CigarUtils.collapseCigarOps;
import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.BASE_MODIFICATIONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_BYTES;
import static com.hartwig.hmftools.common.codon.Nucleotides.baseIndex;
import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;
import static com.hartwig.hmftools.common.sequencing.BiomodalBamUtils.MODC_ANNOTATION;
import static com.hartwig.hmftools.common.sequencing.BiomodalBamUtils.getMMValueFromModCReadIndices;
import static com.hartwig.hmftools.common.sequencing.BiomodalBamUtils.getModCReadIndices;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.DUPLEX_ERROR_QUAL;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.SIMPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SequencingType.BIOMODAL;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;
import static com.hartwig.hmftools.redux.consensus.BaseBuilder.INVALID_POSITION;
import static com.hartwig.hmftools.redux.consensus.BaseBuilder.NO_BASE;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.ALIGNMENT_ONLY;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_MATCH;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_MISMATCH;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.H;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.redux.BaseQualAdjustment;
import com.hartwig.hmftools.common.sequencing.SequencingType;

import org.apache.commons.lang3.Validate;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;
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

        if(sequencingType == BIOMODAL)
            return new BiomodalBuilder(refGenome);

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

        @Override
        public String toString()
        {
            return RefPos + "." + InsertIndex;
        }
    }

    @VisibleForTesting
    public static class AnnotatedBase
    {
        public final ExtendedRefPos Pos;
        public final byte Base;
        public final byte Qual;
        public final Set<String> Annotations;

        private CigarOperator mCigarOp;

        public AnnotatedBase(final ExtendedRefPos pos, byte base, byte qual, final CigarOperator cigarOp)
        {
            Pos = pos;
            Base = base;
            Qual = qual;
            mCigarOp = cigarOp;
            Annotations = Sets.newHashSet();
        }

        public CigarOperator cigarOp() { return mCigarOp; }
        public void softclip()
        {
            Validate.isTrue(mCigarOp.consumesReadBases());
            mCigarOp = S;
        }

        @Override
        public String toString()
        {
            return "AnnotatedBase{" +
                    "Pos=" + Pos +
                    ", Base=" + (char) Base +
                    ", Qual=" + Qual +
                    ", Annotations=" + Annotations +
                    ", CigarOp=" + mCigarOp +
                    '}';
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
                CigarOperator op = el.getOperator() == S ? S : M;
                insertIndex = 0;
                for(int i = 0; i < el.getLength(); i++)
                {
                    ExtendedRefPos pos = new ExtendedRefPos(refPos, insertIndex);
                    annotatedBases.add(new AnnotatedBase(pos, bases[readIndex], quals[readIndex], op));

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
                    annotatedBases.add(new AnnotatedBase(pos, bases[readIndex], quals[readIndex], I));
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

    private static List<List<AnnotatedBase>> annotateReads(final List<SAMRecord> reads)
    {
        return reads.stream().map(NonStandardBaseBuilder::getAnnotatedBases).collect(Collectors.toList());
    }

    @VisibleForTesting
    public static Collection<List<AnnotatedBase>> alignAnnotatedReads(final List<List<AnnotatedBase>> annotateReads)
    {
        SortedMap<ExtendedRefPos, List<AnnotatedBase>> alignment = Maps.newTreeMap();
        for(List<AnnotatedBase> read : annotateReads)
        {
            for(AnnotatedBase base : read)
                alignment.computeIfAbsent(base.Pos, key -> Lists.newArrayList());
        }

        List<ExtendedRefPos> positions = Lists.newArrayList(alignment.keySet());
        for(List<AnnotatedBase> read : annotateReads)
        {
            ExtendedRefPos readPos = read.get(0).Pos;
            int posIdx = 0;
            while(!positions.get(posIdx).equals(readPos))
                posIdx++;

            ExtendedRefPos pos = positions.get(posIdx++);
            alignment.get(pos).add(read.get(0));
            for(int readIdx = 1; readIdx < read.size(); readIdx++)
            {
                readPos = read.get(readIdx).Pos;
                while(!positions.get(posIdx).equals(readPos))
                {
                    pos = positions.get(posIdx++);
                    if(pos.InsertIndex == 0)
                    {
                        alignment.get(pos).add(new AnnotatedBase(pos, NO_BASE, (byte) 0, D));
                        continue;
                    }

                    alignment.get(pos).add(new AnnotatedBase(pos, NO_BASE, (byte) 0, I));
                }

                pos = positions.get(posIdx++);
                alignment.get(pos).add(read.get(readIdx));
            }
        }

        return alignment.values();
    }

    private static List<AnnotatedBase> getConsensusBases(
            final Collection<List<AnnotatedBase>> alignment, final Function<List<AnnotatedBase>, AnnotatedBase> determineConsensus)
    {
        final List<AnnotatedBase> consensus = Lists.newArrayList();
        Consumer<AnnotatedBase> addBase = base ->
        {
            if(base != null)
                consensus.add(base);
        };

        for(List<AnnotatedBase> alignedBases : alignment)
        {
            ExtendedRefPos pos = alignedBases.get(0).Pos;
            if(pos.InsertIndex > 0)
            {
                List<AnnotatedBase> nonMissingBases = Lists.newArrayList();
                int noBaseCount = 0;
                for(AnnotatedBase base : alignedBases)
                {
                    if(base.Base == NO_BASE)
                        noBaseCount++;
                    else
                        nonMissingBases.add(base);
                }

                if(2 * noBaseCount >= alignedBases.size())
                    continue;

                addBase.accept(determineConsensus.apply(nonMissingBases));
                continue;
            }

            List<AnnotatedBase> nonDelBases = Lists.newArrayList();
            int delCount = 0;
            for(AnnotatedBase base : alignedBases)
            {
                if(base.cigarOp() == D)
                {
                    delCount++;
                    continue;
                }

                nonDelBases.add(base);
            }

            if(2 * delCount > alignedBases.size())
                continue;

            addBase.accept(determineConsensus.apply(nonDelBases));
        }

        return consensus;
    }

    private static CigarOperator getConsensusCigarOp(final List<AnnotatedBase> bases)
    {
        if(bases.get(0).Pos.InsertIndex > 0)
            return I;

        for(AnnotatedBase base : bases)
        {
            CigarOperator op = base.cigarOp();
            if(op == M)
                return M;
        }

        return S;
    }

    private static List<AnnotatedBase> insertDels(final List<AnnotatedBase> consensusBases)
    {
        List<AnnotatedBase> newConsensusBases = Lists.newArrayList();
        int lastRefPos = INVALID_POSITION;
        for(AnnotatedBase base : consensusBases)
        {
            int refPos = base.Pos.RefPos;
            if(lastRefPos != INVALID_POSITION && refPos > lastRefPos + 1)
            {
                for(int i = lastRefPos + 1; i < refPos; i++)
                    newConsensusBases.add(new AnnotatedBase(base.Pos, NO_BASE, (byte) 0, D));
            }

            newConsensusBases.add(base);
            lastRefPos = refPos;
        }

        return newConsensusBases;
    }

    private static List<AnnotatedBase> optimizeSoftClips(final RefGenome refGenome, final String chromosome,
            final List<AnnotatedBase> consensusBases)
    {
        List<AnnotatedBase> leftSoftclipBases = Lists.newArrayList();
        List<AnnotatedBase> rightSoftclipBases = Lists.newArrayList();
        List<List<AnnotatedBase>> alignedClusters = Lists.newArrayList();
        List<AnnotatedBase> lastAlignedCluster = null;
        for(AnnotatedBase base : consensusBases)
        {
            if(base.cigarOp() == S)
            {
                if(lastAlignedCluster == null)
                    leftSoftclipBases.add(base);
                else
                    rightSoftclipBases.add(base);

                continue;
            }

            if(lastAlignedCluster == null)
            {
                lastAlignedCluster = Lists.newArrayList(base);
                alignedClusters.add(lastAlignedCluster);
                continue;
            }

            CigarOperator lastCigarOp = lastAlignedCluster.get(0).cigarOp();
            if((lastCigarOp == I || lastCigarOp == D) && lastCigarOp == base.cigarOp())
            {
                lastAlignedCluster.add(base);
                continue;
            }

            lastAlignedCluster = Lists.newArrayList(base);
            alignedClusters.add(lastAlignedCluster);
        }

        if(alignedClusters.isEmpty())
            return consensusBases;

        int[] cScore = new int[alignedClusters.size() + 1];
        cScore[0] = 0;
        for(int i = 0; i < alignedClusters.size(); i++)
        {
            List<AnnotatedBase> cluster = alignedClusters.get(i);
            AnnotatedBase firstBase = cluster.get(0);
            CigarOperator op = firstBase.cigarOp();
            if(op == I || op == D)
            {
                int gapScore = -BWA_GAP_OPEN_PENALTY - cluster.size() * BWA_GAP_EXTEND_PENALTY;
                cScore[i + 1] = cScore[i] + gapScore;
                continue;
            }

            byte refBase = refGenome.getRefBase(chromosome, firstBase.Pos.RefPos);
            int score = refBase == firstBase.Base ? BWA_MATCH_SCORE : -BWA_MISMATCH_PENALTY;
            cScore[i + 1] = cScore[i] + score;
        }

        int minCScore = 0;
        int minCIdx = 0;
        int maxStart = 0;
        int maxEnd = 1;
        int maxScore = cScore[1];
        for(int i = 1; i < cScore.length; i++)
        {
            if(i > 1)
            {
                int score = cScore[i] - minCScore;
                if(score >= maxScore)
                {
                    maxScore = score;
                    maxStart = minCIdx;
                    maxEnd = i;
                }
            }

            if(cScore[i] < minCScore)
            {
                minCScore = cScore[i];
                minCIdx = i;
            }
        }

        boolean extendLeftSoftclip = !leftSoftclipBases.isEmpty();
        if(!extendLeftSoftclip)
        {
            int softClipSaving = -cScore[maxStart];
            if(softClipSaving > BWA_CLIPPING_PENALTY)
                extendLeftSoftclip = true;
        }

        boolean extendRightSoftclip = !rightSoftclipBases.isEmpty();
        if(!extendRightSoftclip)
        {
            int softClipSaving = -(cScore[cScore.length - 1] - cScore[maxEnd]);
            if(softClipSaving > BWA_CLIPPING_PENALTY)
                extendRightSoftclip = true;
        }

        List<AnnotatedBase> finalBases = leftSoftclipBases;
        for(int i = 0; i < alignedClusters.size(); i++)
        {
            List<AnnotatedBase> cluster = alignedClusters.get(i);
            if((i < maxStart && extendLeftSoftclip) || (i >= maxEnd && extendRightSoftclip))
            {
                if(cluster.get(0).cigarOp() == D)
                    continue;

                cluster.forEach(AnnotatedBase::softclip);
            }

            finalBases.addAll(cluster);
        }

        finalBases.addAll(rightSoftclipBases);
        return finalBases;
    }

    private static void finalizeConsensusState(final RefGenome refGenome, final List<SAMRecord> reads, final ConsensusState consensusState,
            boolean hasIndels, final List<AnnotatedBase> consensusBases)
    {
        String chromosome = reads.get(0).getReferenceName();
        final int alignmentStartBoundary = reads.stream().mapToInt(SAMRecord::getAlignmentStart).min().orElse(INVALID_POSITION);
        Validate.isTrue(alignmentStartBoundary >= 1);

        List<AnnotatedBase> finalBases = insertDels(consensusBases);
        finalBases = optimizeSoftClips(refGenome, chromosome, finalBases);

        int minAlignedPosStart = INVALID_POSITION;
        int maxAlignedPosEnd = 0;

        List<Byte> bases = new ArrayList<>(finalBases.size());
        List<Byte> quals = new ArrayList<>(finalBases.size());
        List<CigarOperator> cigarOps = new ArrayList<>(finalBases.size());
        int firstNonClippedIndex = INVALID_POSITION;
        int lastNonClippedIndex = INVALID_POSITION;
        int cigarIdx = 0;
        for(AnnotatedBase base : finalBases)
        {
            int refPos = base.Pos.RefPos;
            CigarOperator op = base.cigarOp();
            cigarOps.add(op);
            if(op == D)
            {
                cigarIdx++;
                continue;
            }

            bases.add(base.Base);
            quals.add(base.Qual);
            if(firstNonClippedIndex == INVALID_POSITION && !op.isClipping())
                firstNonClippedIndex = cigarIdx;

            if(!op.isClipping())
                lastNonClippedIndex = cigarIdx;

            if(minAlignedPosStart == INVALID_POSITION && op == M)
            {
                Validate.isTrue(refPos >= alignmentStartBoundary);
                minAlignedPosStart = refPos;
            }

            if(op == M)
                maxAlignedPosEnd = refPos;

            cigarIdx++;
        }

        for(int i = firstNonClippedIndex; i <= lastNonClippedIndex; i++)
        {
            if(cigarOps.get(i).isClipping())
                cigarOps.set(i, M);
        }

        Cigar cigar = new Cigar(collapseCigarOps(cigarOps));
        int minUnclippedPosStart = minAlignedPosStart - leftSoftClipLength(cigar);
        int maxUnclippedPosEnd = maxAlignedPosEnd + rightSoftClipLength(cigar);
        consensusState.setBoundaries(minUnclippedPosStart, maxUnclippedPosEnd, minAlignedPosStart, maxAlignedPosEnd);
        consensusState.setBaseLength(bases.size());

        for(int i = 0; i < bases.size(); i++)
        {
            consensusState.Bases[i] = bases.get(i);
            consensusState.BaseQualities[i] = quals.get(i);
        }

        consensusState.CigarElements.addAll(cigar.getCigarElements());

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

            List<List<AnnotatedBase>> annotatedReads = annotateReads(reads);
            Collection<List<AnnotatedBase>> alignment = alignAnnotatedReads(annotatedReads);
            List<AnnotatedBase> consensusBases = getConsensusBases(alignment, records -> determineConsensus(chromosome, records));
            finalizeConsensusState(mRefGenome, reads, consensusState, hasIndels, consensusBases);
        }

        @VisibleForTesting
        @Nullable
        public AnnotatedBase determineConsensus(final String chromosome, final List<AnnotatedBase> bases)
        {
            ExtendedRefPos pos = bases.get(0).Pos;
            int refPos = pos.RefPos;
            CigarOperator consensusOp = getConsensusCigarOp(bases);
            if(consensusOp != M || refPos < 1 || refPos > mChromosomeLength)
                refPos = INVALID_POSITION;

            if(bases.size() == 1)
            {
                boolean isQualZero = bases.get(0).Qual == (byte) 0;
                if(isQualZero && refPos != INVALID_POSITION)
                {
                    byte refBase = mRefGenome.getRefBase(chromosome, refPos);
                    return new AnnotatedBase(pos, refBase, (byte) 1, M);
                }

                if(isQualZero && consensusOp == I)
                    return null;

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
                if(refPos != INVALID_POSITION)
                {
                    byte refBase = mRefGenome.getRefBase(chromosome, refPos);
                    return new AnnotatedBase(pos, refBase, (byte) 1, M);
                }

                if(consensusOp == I)
                    return null;

                return new AnnotatedBase(pos, ANY_BASE, (byte) 0, S);
            }

            if(!baseCountsByQual.containsKey((byte) SIMPLEX_QUAL) && !baseCountsByQual.containsKey((byte) DUPLEX_QUAL))
            {
                if(refPos != INVALID_POSITION)
                {
                    byte refBase = mRefGenome.getRefBase(chromosome, refPos);
                    return new AnnotatedBase(pos, refBase, (byte) 1, M);
                }

                if(consensusOp == I)
                    return null;

                return new AnnotatedBase(pos, ANY_BASE, (byte) 0, S);
            }

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

                if(refPos != INVALID_POSITION)
                {
                    byte refBase = mRefGenome.getRefBase(chromosome, refPos);
                    return new AnnotatedBase(pos, refBase, (byte) 1, M);
                }

                if(consensusOp == I)
                    return null;

                return new AnnotatedBase(pos, ANY_BASE, (byte) 0, S);
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
                    multipleMax = true;
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
                if(refPos != INVALID_POSITION)
                {
                    byte refBase = mRefGenome.getRefBase(chromosome, refPos);
                    return new AnnotatedBase(pos, refBase, (byte) 1, M);
                }

                if(consensusOp == I)
                    return null;

                return new AnnotatedBase(pos, ANY_BASE, (byte) 0, S);
            }

            byte base = DNA_BASE_BYTES[maxIdx];
            byte qual = (byte) SIMPLEX_QUAL;
            return new AnnotatedBase(pos, base, qual, consensusOp);
        }
    }

    public static class BiomodalBuilder extends NonStandardBaseBuilder
    {
        private final BaseBuilder mBaseBuilder;

        public BiomodalBuilder(final RefGenome refGenome)
        {
            super(refGenome);

            mBaseBuilder = new BaseBuilder(refGenome, null);
        }

        @Override
        public void buildConsensusRead(final List<SAMRecord> reads, final ConsensusState consensusState, final boolean hasIndels)
        {
            String chromosome = reads.get(0).getReferenceName();
            mChromosomeLength = mRefGenome.getChromosomeLength(chromosome);

            List<List<AnnotatedBase>> annotatedReads = annotateReads(reads);
            for(int i = 0; i < reads.size(); i++)
                addModCAnnotation(reads.get(i), annotatedReads.get(i));

            Collection<List<AnnotatedBase>> alignment = alignAnnotatedReads(annotatedReads);
            List<AnnotatedBase> consensusBases = getConsensusBases(alignment, records -> determineConsensus(chromosome, consensusState.IsForward, records));
            finalizeConsensusState(mRefGenome, reads, consensusState, hasIndels, consensusBases);

            SortedSet<Integer> modCReadIndices = Sets.newTreeSet();
            for(int i = 0; i < consensusBases.size(); i++)
            {
                if(consensusBases.get(i).Annotations.contains(MODC_ANNOTATION))
                    modCReadIndices.add(i);
            }

            String mmValue = getMMValueFromModCReadIndices(consensusState.Bases, modCReadIndices, consensusState.IsForward);
            consensusState.Attributes.put(BASE_MODIFICATIONS_ATTRIBUTE, mmValue);
        }

        private AnnotatedBase determineConsensus(final String chromosome, boolean isForward, final List<AnnotatedBase> bases)
        {
            CigarOperator consensusOp = getConsensusCigarOp(bases);
            ExtendedRefPos pos = bases.get(0).Pos;
            byte[] locationBases = new byte[bases.size()];
            byte[] locationQuals = new byte[bases.size()];
            for(int i = 0; i < bases.size(); i++)
            {
                locationBases[i] = bases.get(i).Base;
                locationQuals[i] = bases.get(i).Qual;
            }

            int basePosition = pos.RefPos;
            if(consensusOp != M || basePosition < 1 || basePosition > mChromosomeLength)
                basePosition = INVALID_POSITION;

            byte[] consensusBaseAndQual = mBaseBuilder.determineBaseAndQual(locationBases, locationQuals, chromosome, basePosition);
            byte consensusBase = consensusBaseAndQual[0] == NO_BASE ? ANY_BASE : consensusBaseAndQual[0];
            byte consensusQual = BaseQualAdjustment.adjustBaseQual(consensusBaseAndQual[1]);

            if(consensusQual == (byte) 0 && consensusOp == I)
                return null;

            AnnotatedBase consensus = new AnnotatedBase(pos, consensusBase, consensusQual, consensusOp);

            byte targetBase = isForward ? (byte) 'C' : (byte) swapDnaBase('C');
            if(consensusBase != targetBase)
                return consensus;

            int nonModCQualTotal = 0;
            int modCQualTotal = 0;
            for(AnnotatedBase annotatedBase : bases)
            {
                if(annotatedBase.Base != targetBase)
                    continue;

                if(annotatedBase.Annotations.contains(MODC_ANNOTATION))
                {
                    modCQualTotal += annotatedBase.Qual;
                    continue;
                }

                nonModCQualTotal += annotatedBase.Qual;
            }

            if(modCQualTotal <= nonModCQualTotal)
                return consensus;

            consensus.Annotations.add(MODC_ANNOTATION);
            return consensus;
        }

        private static void addModCAnnotation(final SAMRecord read, List<AnnotatedBase> annotatedRead)
        {
            SortedSet<Integer> modCReadIndices = getModCReadIndices(read);
            for(int readIndex : modCReadIndices)
                annotatedRead.get(readIndex).Annotations.add(MODC_ANNOTATION);
        }
    }
}
