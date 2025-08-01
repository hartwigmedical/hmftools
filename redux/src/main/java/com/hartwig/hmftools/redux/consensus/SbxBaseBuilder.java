package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_BYTES;
import static com.hartwig.hmftools.common.codon.Nucleotides.baseIndex;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SIMPLEX_QUAL;
import static com.hartwig.hmftools.redux.consensus.BaseBuilder.INVALID_POSITION;

import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.Collection;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class SbxBaseBuilder extends NonStandardBaseBuilder
{
    public static final int DUPLEX_ERROR_QUAL = 0;

    public SbxBaseBuilder(final RefGenome refGenome)
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
        {
            refPos = INVALID_POSITION;
        }

        if(bases.size() == 1)
        {
            boolean isQualZero = bases.get(0).Qual == (byte) 0;
            if(isQualZero && refPos != INVALID_POSITION)
            {
                byte refBase = mRefGenome.getRefBase(chromosome, refPos);
                return new AnnotatedBase(pos, refBase, (byte) 1, M);
            }

            if(isQualZero && consensusOp == I)
            {
                return null;
            }

            return bases.get(0);
        }

        Map<Byte, int[]> baseCountsByQual = Maps.newHashMap();
        for(AnnotatedBase annotatedBase : bases)
        {
            byte base = annotatedBase.Base;
            int baseIdx = baseIndex(base);
            if(baseIdx < 0 || baseIdx >= DNA_BASE_BYTES.length)
            {
                continue;
            }

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
            {
                return null;
            }

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
            {
                return null;
            }

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
            {
                return null;
            }

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
            {
                multipleMax = true;
            }
        }

        int[] duplexErrorCounts = baseCountsByQual.get((byte) DUPLEX_ERROR_QUAL);
        if(duplexErrorCounts != null)
        {
            for(int duplexErrorCount : duplexErrorCounts)
            {
                totalCount += duplexErrorCount;
            }
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
            {
                return null;
            }

            return new AnnotatedBase(pos, ANY_BASE, (byte) 0, S);
        }

        byte base = DNA_BASE_BYTES[maxIdx];
        byte qual = (byte) SIMPLEX_QUAL;
        return new AnnotatedBase(pos, base, qual, consensusOp);
    }
}
