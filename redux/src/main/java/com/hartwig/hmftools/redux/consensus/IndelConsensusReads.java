package com.hartwig.hmftools.redux.consensus;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.maxQual;
import static com.hartwig.hmftools.redux.ReduxConfig.isSbx;
import static com.hartwig.hmftools.redux.consensus.SbxIndelConsensus.determineIndelConsensus;
import static com.hartwig.hmftools.redux.duplicate.DuplicateGroupBuilder.calcBaseQualAverage;
import static com.hartwig.hmftools.redux.consensus.BaseQualPair.NO_BASE;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_FAIL;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_MATCH;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_MISMATCH;
import static com.hartwig.hmftools.redux.consensus.ConsensusState.alignedOrClipped;
import static com.hartwig.hmftools.redux.consensus.ConsensusState.consumesRefOrUnclippedBases;
import static com.hartwig.hmftools.redux.consensus.ConsensusState.deleteOrSplit;
import static com.hartwig.hmftools.redux.consensus.IlluminaRoutines.isDualStrandAndIsFirstInPair;

import static htsjdk.samtools.CigarOperator.H;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.redux.BaseQualAdjustment;

import htsjdk.samtools.CigarElement;
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
            consensusState.setFromRead(templateRead, false);
            mBaseBuilder.buildReadBases(reads, consensusState);
            consensusState.setOutcome(INDEL_MATCH);
            return;
        }

        if(isSbx())
        {
            determineIndelConsensus(this, consensusState, reads);
            return;
        }

        int readCount = reads.size();

        boolean[] isFirstInPair = new boolean[readCount];
        boolean isDualStrand = isDualStrandAndIsFirstInPair(reads, isFirstInPair);

        // find the most common read by CIGAR, and where there are equal counts choose the one with the least soft-clips
        int baseLength = templateRead.getReadBases().length;
        consensusState.setBaseLength(baseLength);
        consensusState.setBoundaries(templateRead);

        List<ReadParseState> readStates = reads.stream().map(x -> new ReadParseState(x, consensusState.IsForward)).collect(Collectors.toList());

        int baseIndex = consensusState.IsForward ? 0 : baseLength - 1;
        int refPosition = consensusState.IsForward ? consensusState.UnclippedPosStart : consensusState.UnclippedPosEnd;

        int initialRefPosition = refPosition;
        readStates.forEach(x -> x.moveToRefPosition(initialRefPosition));

        List<CigarElement> selectElements = templateRead.getCigar().getCigarElements();
        int cigarCount = selectElements.size();
        int cigarIndex = consensusState.IsForward ? 0 : cigarCount - 1;

        while(cigarIndex >= 0 && cigarIndex < cigarCount)
        {
            CigarElement element = selectElements.get(cigarIndex);

            // simplest scenario is where all reads agree about this next element
            addElementBases(consensusState, readStates, element, baseIndex, refPosition, isDualStrand, isFirstInPair);

            if(consensusState.outcome() == INDEL_FAIL)
                break;

            if(element.getOperator().consumesReadBases())
            {
                if(consensusState.IsForward)
                    baseIndex += element.getLength();
                else
                    baseIndex -= element.getLength();
            }

            if(consumesRefOrUnclippedBases(element.getOperator()))
            {
                if(consensusState.IsForward)
                    refPosition += element.getLength();
                else
                    refPosition -= element.getLength();
            }

            if(consensusState.IsForward)
                ++cigarIndex;
            else
                --cigarIndex;
        }

        consensusState.finaliseCigar();

        if(consensusState.outcome() != INDEL_FAIL)
            consensusState.setOutcome(INDEL_MISMATCH);
    }

    protected void addElementBases(
            final ConsensusState consensusState, final List<ReadParseState> readStates, final CigarElement selectedElement,
            int baseIndex, int refPosition, boolean isDualStrand, final boolean[] isFirstInPair)
    {
        int chromosomeLength = mBaseBuilder.chromosomeLength();
        if(chromosomeLength == 0)
            chromosomeLength = mBaseBuilder.refGenome().getChromosomeLength(readStates.get(0).Read.getReferenceName());

        int readCount = readStates.size();

        consensusState.addCigarElement(selectedElement.getLength(), selectedElement.getOperator());

        if(deleteOrSplit(selectedElement.getOperator()) || selectedElement.getOperator() == H)
        {
            int newRefPosition = refPosition + (consensusState.IsForward ? 1 : -1) * selectedElement.getLength();

            // move past the delete element and any differing aligned bases
            for(int r = 0; r < readCount; ++r)
            {
                ReadParseState read = readStates.get(r);

                if(read.exhausted() || read.beforeUnclippedPosition(newRefPosition))
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

                    if(deleteOrSplit(read.elementType()))
                    {
                        read.moveNextBase();
                    }
                    else if(alignedOrClipped(read.elementType()))
                    {
                        if(read.refPosition() == newRefPosition)
                            break;

                        read.moveNextBase();
                    }
                }
            }

            return;
        }

        byte[] locationBases = new byte[readCount];
        byte[] locationQuals = new byte[readCount];
        int newRefPosition = refPosition;

        for(int i = 0; i < selectedElement.getLength(); ++i)
        {
            boolean hasMismatch = false;
            byte maxQual = BaseQualAdjustment.BASE_QUAL_MINIMUM;
            byte firstBase = NO_BASE;

            if(i > 0) // reset consensus arrays
            {
                for(int r = 0; r < readCount; ++r)
                {
                    locationBases[r] = NO_BASE;
                    locationQuals[r] = 0;
                }
            }

            for(int r = 0; r < readCount; ++r)
            {
                ReadParseState read = readStates.get(r);

                if(read.exhausted() || read.beforeUnclippedPosition(newRefPosition))
                    continue;

                // check for element type differences:

                // first skip past any insert if the selected element is aligned
                if(selectedElement.getOperator().consumesReferenceBases() && read.elementType() == I)
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
                        if(read.elementType() == M || deleteOrSplit(read.elementType()) || read.elementType() == S)
                        {
                            moveNext = false;
                            useBase = false;
                        }
                    }
                }
                else if(selectedElement.getOperator() == I && read.elementLength() != selectedElement.getLength() && isSbx())
                {
                    // only apply this additional check to SBX to force it to use quals of reads matching the selected insert exactly
                    useBase = false;
                }

                if(useBase)
                {
                    locationBases[r] = read.base();
                    locationQuals[r] = read.baseQual();

                    if(firstBase == NO_BASE)
                        firstBase = locationBases[r];
                    else
                        hasMismatch |= locationBases[r] != firstBase;

                    maxQual = maxQual(locationQuals[r], maxQual);
                }

                if(moveNext)
                    read.moveNextBase();
            }

            if(!hasMismatch && firstBase != NO_BASE)
            {
                consensusState.Bases[baseIndex] = firstBase;
                consensusState.BaseQualities[baseIndex] = maxQual;
            }
            else
            {
                int basePosition = newRefPosition < 1 || newRefPosition > chromosomeLength ? BaseBuilder.INVALID_POSITION : newRefPosition;

                BaseQualPair consensusBaseAndQual = mBaseBuilder.determineBaseAndQual(
                        locationBases, locationQuals, consensusState.Chromosome, basePosition, isDualStrand, isFirstInPair);

                consensusState.Bases[baseIndex] = consensusBaseAndQual.Base;
                consensusState.BaseQualities[baseIndex] = consensusBaseAndQual.Qual;
            }

            if(consensusState.IsForward)
            {
                ++baseIndex;

                if(consumesRefOrUnclippedBases(selectedElement.getOperator()))
                    ++newRefPosition;
            }
            else
            {
                --baseIndex;

                if(consumesRefOrUnclippedBases(selectedElement.getOperator()))
                    --newRefPosition;
            }
        }
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
