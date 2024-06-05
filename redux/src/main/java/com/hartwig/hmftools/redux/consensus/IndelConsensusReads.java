package com.hartwig.hmftools.redux.consensus;

import static java.lang.Math.max;

import static com.hartwig.hmftools.redux.common.DuplicateGroupBuilder.calcBaseQualAverage;
import static com.hartwig.hmftools.redux.consensus.BaseBuilder.INVALID_POSITION;
import static com.hartwig.hmftools.redux.consensus.BaseBuilder.NO_BASE;
import static com.hartwig.hmftools.redux.consensus.BaseBuilder.isDualStrandAndIsFirstInPair;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_FAIL;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_MATCH;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_MISMATCH;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.H;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.N;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.qual.BaseQualAdjustment;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class IndelConsensusReads
{
    private final BaseBuilder mBaseBuilder;

    public IndelConsensusReads(final BaseBuilder baseBuilder)
    {
        mBaseBuilder = baseBuilder;
    }

    public void buildIndelComponents(final List<SAMRecord> reads, final ConsensusState consensusState, final SAMRecord templateRead)
    {
        boolean hasCigarMismatch = reads.stream().anyMatch(x -> !x.getCigarString().equals(templateRead.getCigarString()));

        if(!hasCigarMismatch)
        {
            // SAMRecord selectedConsensusRead = reads.get(0);
            int baseLength = templateRead.getReadBases().length;
            consensusState.setBaseLength(baseLength);
            consensusState.setBoundaries(templateRead);

            mBaseBuilder.buildReadBases(reads, consensusState);
            consensusState.setOutcome(INDEL_MATCH);
            consensusState.CigarElements.addAll(templateRead.getCigar().getCigarElements());
            return;
        }

        int readCount = reads.size();

        boolean[] isFirstInPair = new boolean[readCount];
        boolean isDualStrand = isDualStrandAndIsFirstInPair(reads, isFirstInPair);

        // find the most common read by CIGAR, and where there are equal counts choose the one with the least soft-clips
        // SAMRecord selectedConsensusRead = selectConsensusRead(cigarFrequencies);
        // SAMRecord selectedConsensusRead = templateRead;

        int baseLength = templateRead.getReadBases().length;
        consensusState.setBaseLength(baseLength);
        consensusState.setBoundaries(templateRead);

        List<ReadParseState> readStates = reads.stream().map(x -> new ReadParseState(x, consensusState.IsForward)).collect(Collectors.toList());

        int baseIndex = consensusState.IsForward ? 0 : baseLength - 1;

        List<CigarElement> selectElements = templateRead.getCigar().getCigarElements();
        int cigarCount = selectElements.size();
        int cigarIndex = consensusState.IsForward ? 0 : cigarCount - 1;

        while(cigarIndex >= 0 && cigarIndex < cigarCount)
        {
            CigarElement element = selectElements.get(cigarIndex);

            // simplest scenario is where all reads agree about this next element
            addElementBases(consensusState, readStates, element, baseIndex, isDualStrand, isFirstInPair);

            if(consensusState.outcome() == INDEL_FAIL)
                break;

            if(element.getOperator().consumesReadBases())
            {
                if(consensusState.IsForward)
                    baseIndex += element.getLength();
                else
                    baseIndex -= element.getLength();
            }

            if(consensusState.IsForward)
                ++cigarIndex;
            else
                --cigarIndex;
        }

        if(consensusState.outcome() != INDEL_FAIL)
            consensusState.setOutcome(INDEL_MISMATCH);
    }

    private void addElementBases(
            final ConsensusState consensusState, final List<ReadParseState> readStates, final CigarElement selectedElement, int baseIndex,
            boolean isDualStrand, final boolean[] isFirstInPair)
    {
        int chromosomeLength = mBaseBuilder.chromosomeLength();
        if(chromosomeLength == 0)
            chromosomeLength = mBaseBuilder.refGenome().getChromosomeLength(readStates.get(0).Read.getReferenceName());

        int readCount = readStates.size();

        consensusState.addCigarElement(selectedElement.getLength(), selectedElement.getOperator());

        if(deleteOrSplit(selectedElement.getOperator()) || selectedElement.getOperator() == H)
        {
            // move past the delete element and any differing aligned bases
            for(int r = 0; r < readCount; ++r)
            {
                ReadParseState read = readStates.get(r);

                if(read.exhausted())
                    continue;

                if(read.elementType() == selectedElement.getOperator()
                && read.elementIndex() == 0 && read.elementLength() == selectedElement.getLength())
                {
                    read.skipNonReadBaseElement();
                    continue;
                }

                for(int i = 0; i < selectedElement.getLength(); ++i)
                {
                    if(read.elementType() == I)
                        read.skipInsertElement();

                    if(deleteOrSplit(read.elementType()) || alignedOrClipped(read.elementType()))
                    {
                        read.moveNextBase();
                    }
                }
            }

            return;
        }

        byte[] locationBases = new byte[readCount];
        byte[] locationQuals = new byte[readCount];

        for(int i = 0; i < selectedElement.getLength(); ++i)
        {
            boolean hasMismatch = false;
            int maxQual = 0;
            byte firstBase = NO_BASE;

            for(int r = 0; r < readCount; ++r)
            {
                locationBases[r] = NO_BASE;
            }

            for(int r = 0; r < readCount; ++r)
            {
                ReadParseState read = readStates.get(r);

                if(read.exhausted())
                    continue;

                // check for element type differences:

                // first skip past any insert if the selected element is aligned
                if(selectedElement.getOperator() == M && read.elementType() == I)
                    read.skipInsertElement();

                boolean useBase = true;
                boolean moveNext = true;

                if(read.elementType() != selectedElement.getOperator())
                {
                    // when aligned (M) is selected:
                    // - insert - skip past the insert's bases
                    // - delete - move along in step but cannot use base
                    // - ignore read with soft-clipped bases since they are likely misaligned

                    // when insert (I) is selected:
                    // - aligned or delete - pause the index and don't use the base

                    if(selectedElement.getOperator() == M)
                    {
                        if(deleteOrSplit(read.elementType()) || read.elementType() == S)
                        {
                            useBase = false;
                        }
                        else if(read.elementType() == I)
                        {
                            // handled above, implies a bug or consecutive insert
                            // logMismatchFail(consensusState, r, read, selectedElement);
                            consensusState.setOutcome(INDEL_FAIL);
                            return;
                        }
                    }
                    else if(selectedElement.getOperator() == I)
                    {
                        if(read.elementType() == M || deleteOrSplit(read.elementType()))
                        {
                            moveNext = false;
                            useBase = false;
                        }
                    }
                }

                if(useBase)
                {
                    locationBases[r] = read.currentBase();
                    locationQuals[r] = read.currentBaseQual();

                    if(firstBase == NO_BASE)
                        firstBase = locationBases[r];
                    else
                        hasMismatch |= locationBases[r] != firstBase;

                    maxQual = max(locationQuals[r], maxQual);
                }

                if(moveNext)
                    read.moveNextBase();
            }

            if(!hasMismatch)
            {
                consensusState.Bases[baseIndex] = firstBase;
                consensusState.BaseQualities[baseIndex] = (byte)maxQual;
            }
            else
            {
                int basePosition = consensusState.MinUnclippedPosStart + baseIndex;

                if(basePosition < 1 || basePosition > chromosomeLength)
                    basePosition = BaseBuilder.INVALID_POSITION;

                byte[] consensusBaseAndQual;

                if(isDualStrand && basePosition != INVALID_POSITION)
                {
                    // split the reads into 2 consensus reads and then compare
                    consensusBaseAndQual = mBaseBuilder.determineDualStrandBaseAndQual(
                            isFirstInPair, locationBases, locationQuals, consensusState.Chromosome, basePosition);
                }
                else
                {
                    consensusBaseAndQual = mBaseBuilder.determineBaseAndQual(
                            locationBases, locationQuals, consensusState.Chromosome, basePosition);
                }

                consensusState.Bases[baseIndex] = consensusBaseAndQual[0];
                consensusState.BaseQualities[baseIndex] = BaseQualAdjustment.adjustBaseQual(consensusBaseAndQual[1]);
            }

            if(consensusState.IsForward)
                ++baseIndex;
            else
                --baseIndex;
        }
    }

    private static boolean deleteOrSplit(final CigarOperator operator)
    {
        return operator == D || operator == N;
    }

    public static boolean alignedOrClipped(final CigarOperator operator)
    {
        return operator == M || operator.isClipping();
    }

    public static SAMRecord selectPrimaryRead(final List<SAMRecord> reads)
    {
        // choose the read with the longest aligned bases
        int maxAlignedBases = 0;
        double maxBaseQual = 0;
        SAMRecord maxRead = null;

        for(SAMRecord read : reads)
        {
            int alignedBases = read.getCigar().getCigarElements().stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();

            if(alignedBases > maxAlignedBases)
            {
                maxRead = read;
                maxAlignedBases = alignedBases;
                maxBaseQual = calcBaseQualAverage(read);
            }
            else if(alignedBases == maxAlignedBases)
            {
                double avgBaseQual = calcBaseQualAverage(read);

                if(avgBaseQual > maxBaseQual)
                {
                    maxRead = read;
                    maxAlignedBases = alignedBases;
                    maxBaseQual = avgBaseQual;
                }
            }
        }

        return maxRead;
    }
}
