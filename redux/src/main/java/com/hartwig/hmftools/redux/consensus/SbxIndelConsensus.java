package com.hartwig.hmftools.redux.consensus;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.RAW_SIMPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_MISMATCH_QUAL;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_FAIL;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_MISMATCH;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_SOFTCLIP;
import static com.hartwig.hmftools.redux.consensus.ConsensusState.consumesRefOrUnclippedBases;
import static com.hartwig.hmftools.redux.consensus.ConsensusState.deleteOrSplit;

import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Arrays;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class SbxIndelConsensus
{
    public static SAMRecord checkSoftClipIndelPair(final List<SAMRecord> reads)
    {
        if(reads.size() != 2)
            return null;

        // check for a soft-clip overlapping indels and if found use the read with the soft-clips
        for(int i = 0; i <= 1; ++i)
        {
            SAMRecord softClipRead = reads.get(i);
            SAMRecord indelRead = reads.get(1 - i);

            for(int j = 0; j <= 1; ++j)
            {
                boolean checkLeft = (j == 0);

                int softClipLength = checkLeft ? leftSoftClipLength(softClipRead) : rightSoftClipLength(softClipRead);

                if(softClipLength == 0)
                    continue;

                if(softClipOverlapsIndels(softClipLength, indelRead, checkLeft))
                    return softClipRead;
            }
        }

        return null;
    }

    private static boolean softClipOverlapsIndels(int softClipLength, final SAMRecord otherRead, boolean checkLeft)
    {
        int readBases = 0;
        List<CigarElement> cigarElements = otherRead.getCigar().getCigarElements();
        int cigarCount = cigarElements.size();
        int cigarIndex = checkLeft ? 0 : cigarCount - 1;
        while(cigarIndex >= 0 && cigarIndex < cigarCount)
        {
            CigarElement element = cigarElements.get(cigarIndex);

            if(element.getOperator().consumesReadBases())
            {
                if(element.getOperator() == I && readBases <= softClipLength)
                    return true;

                readBases += element.getLength();

                if(readBases > softClipLength)
                    return false;
            }

            if(checkLeft)
                ++cigarIndex;
            else
                --cigarIndex;
        }

        return false;
    }

    public static void determineIndelConsensus(
            final IndelConsensusReads indelConsensusReads, final ConsensusState consensusState, final List<SAMRecord> reads)
    {
        SAMRecord softClipOverIndelRead = checkSoftClipIndelPair(reads);

        if(softClipOverIndelRead != null)
        {
            consensusState.setFromRead(softClipOverIndelRead, true);
            consensusState.setOutcome(INDEL_SOFTCLIP);
            return;
        }

        List<ReadParseState> readStates = reads.stream().map(x -> new ReadParseState(x, consensusState.IsForward)).collect(Collectors.toList());

        int outerUnclippedPosition = -1;

        for(ReadParseState read : readStates)
        {
            if(consensusState.IsForward)
            {
                int readStart = read.Read.getAlignmentStart();

                int readUnclippedPosition = readStart - leftSoftClipLength(read.Read);

                if(outerUnclippedPosition < 0 || readUnclippedPosition < outerUnclippedPosition)
                    outerUnclippedPosition = readUnclippedPosition;
            }
            else
            {
                int readEnd = read.Read.getAlignmentEnd();
                int readUnclippedPosition = readEnd + rightSoftClipLength(read.Read);
                outerUnclippedPosition = max(outerUnclippedPosition, readUnclippedPosition);
            }
        }

        int initialRefPosition = outerUnclippedPosition;
        readStates.forEach(x -> x.moveToRefPosition(initialRefPosition));

        List<CigarElement> cigarElements = determineConsensusCigar(readStates, consensusState.IsForward, initialRefPosition);

        int baseLength = cigarElements.stream().filter(x -> x.getOperator().consumesReadBases()).mapToInt(x -> x.getLength()).sum();
        consensusState.setBaseLength(baseLength);

        int baseIndex = consensusState.IsForward ? 0 : baseLength - 1;
        int refPosition = initialRefPosition;

        int cigarCount = cigarElements.size();
        int cigarIndex = consensusState.IsForward ? 0 : cigarCount - 1;
        boolean[] isFirstInPair = null; // unused

        readStates.forEach(x -> x.reset());
        readStates.forEach(x -> x.moveToRefPosition(initialRefPosition));

        while(cigarIndex >= 0 && cigarIndex < cigarCount)
        {
            CigarElement element = cigarElements.get(cigarIndex);

            indelConsensusReads.addElementBases(consensusState, readStates, element, baseIndex, refPosition, false, isFirstInPair);

            if(consensusState.outcome() == INDEL_FAIL)
                break;

            if(consensusState.IsForward)
            {
                if(element.getOperator().consumesReadBases())
                    baseIndex += element.getLength();

                if(consumesRefOrUnclippedBases(element.getOperator()))
                    refPosition += element.getLength();

                ++cigarIndex;
            }
            else
            {
                if(element.getOperator().consumesReadBases())
                    baseIndex -= element.getLength();

                if(consumesRefOrUnclippedBases(element.getOperator()))
                    refPosition -= element.getLength();

                --cigarIndex;
            }
        }

        consensusState.finaliseCigar();

        int consensusAlignmentStart, consensusAlignmentEnd, consensusUnclippedStart, consensusUnclippedEnd;

        int refBaseLength = cigarElements.stream().filter(x -> x.getOperator().consumesReferenceBases()).mapToInt(x -> x.getLength()).sum();

        if(consensusState.IsForward)
        {
            consensusUnclippedStart = outerUnclippedPosition;
            consensusAlignmentStart = outerUnclippedPosition;

            if(cigarElements.get(0).getOperator() == S)
                consensusAlignmentStart += cigarElements.get(0).getLength();

            consensusAlignmentEnd = consensusAlignmentStart + refBaseLength - 1;

            consensusUnclippedEnd = consensusAlignmentEnd;

            if(cigarElements.get(cigarElements.size() - 1).getOperator() == S)
                consensusUnclippedEnd += cigarElements.get(cigarElements.size() - 1).getLength();
        }
        else
        {
            consensusUnclippedEnd = outerUnclippedPosition;
            consensusAlignmentEnd = consensusUnclippedEnd;

            if(cigarElements.get(cigarElements.size() - 1).getOperator() == S)
                consensusAlignmentEnd -= cigarElements.get(cigarElements.size() - 1).getLength();

            consensusAlignmentStart = consensusAlignmentEnd - refBaseLength + 1;

            consensusUnclippedStart = consensusAlignmentStart;

            if(cigarElements.get(0).getOperator() == S)
                consensusUnclippedStart -= cigarElements.get(0).getLength();
        }

        consensusState.setBoundaries(consensusUnclippedStart, consensusUnclippedEnd, consensusAlignmentStart, consensusAlignmentEnd);

        if(consensusState.outcome() != INDEL_FAIL)
            consensusState.setOutcome(INDEL_MISMATCH);
    }

    private static List<CigarElement> determineConsensusCigar(
            final List<ReadParseState> readStates, boolean isForwardConsensus, int initialUnclippedPosition)
    {
        List<CigarElement> cigarElements = Lists.newArrayList();

        /* routine:
        - determine earliest soft-clip start
        - determine earliest aligned start
        - build left soft-clip element using any read covering that unclipped read base up until the earlier aligned start
            - that favours M over S
        - then begin element prioritisation routine:
	        - walk forwards while read cigars agree
	        - at point of difference, take element type based on rules above
	        - this will then use that read's element until exhausted
	        - move all reads to this same point
	        - repeat process
        */

        List<ReadParseState> activeReadStates = Lists.newArrayList(readStates);

        int elementLength = 0;
        CigarOperator elementOperator = null;
        boolean inAlignedBases = false; // to prevent unvalid changes back from aligned to initial S
        boolean inFinalSoftClip = false; // to prevent invalid changes from the final S back to aligned
        int refPosition = initialUnclippedPosition;

        while(true)
        {
            CigarOperator nextOperator = inFinalSoftClip ?
                    S : determineConsensusCigarOperator(activeReadStates, inAlignedBases, elementOperator, refPosition, isForwardConsensus);

            if(elementOperator == null || elementOperator != nextOperator)
            {
                if(elementOperator != null)
                {
                    if(isForwardConsensus)
                        cigarElements.add(new CigarElement(elementLength, elementOperator));
                    else
                        cigarElements.add(0, new CigarElement(elementLength, elementOperator));
                }

                elementLength = 0;
                elementOperator = nextOperator;

                if(!inAlignedBases && elementOperator == M)
                    inAlignedBases = true;
                else if(inAlignedBases && elementOperator == S)
                    inFinalSoftClip = true;
            }

            ++elementLength;

            // move reads ahead depending on whether they agree with the consensus or not
            boolean hasExhausted = false;

            for(ReadParseState read : activeReadStates)
            {
                if(read.beforeUnclippedPosition(refPosition))
                    continue;

                // check for element type differences

                // first skip past any insert if the selected element is aligned
                if(nextOperator.consumesReferenceBases() && read.elementType() == I)
                    read.skipInsertElement();

                boolean moveNext = true;

                if(read.elementType() != nextOperator)
                {
                    // when aligned (M) is selected:
                    // - insert - skip past the insert's bases
                    // - delete - move along in step but cannot use base
                    // - ignore read with soft-clipped bases since they are likely misaligned

                    // when insert (I) is selected:
                    // - aligned or delete - pause the index and don't use the base
                    if(nextOperator == I)
                    {
                        // keep all other operators fixed, including soft-clips
                        moveNext = false;
                    }
                }

                if(moveNext)
                    read.moveNextBase();

                hasExhausted |= read.exhausted();
            }

            if(hasExhausted)
            {
                int index = 0;
                while(index < activeReadStates.size())
                {
                    if(activeReadStates.get(index).exhausted())
                        activeReadStates.remove(index);
                    else
                        ++index;
                }

                if(activeReadStates.isEmpty())
                    break;
            }

            if(consumesRefOrUnclippedBases(nextOperator))
                refPosition += isForwardConsensus ? 1 : -1;
        }

        // add the last
        if(isForwardConsensus)
            cigarElements.add(new CigarElement(elementLength, elementOperator));
        else
            cigarElements.add(0, new CigarElement(elementLength, elementOperator));

        return cigarElements;
    }

    private static CigarOperator determineConsensusCigarOperator(
            final List<ReadParseState> readStates, boolean inAlignedBases, final CigarOperator currentOperator,
            int refPosition, boolean isForwardConsensus)
    {
        CigarOperator operator = null;
        boolean matched = true;

        for(ReadParseState read : readStates)
        {
            if(read.beforeUnclippedPosition(refPosition))
                continue;

            CigarOperator readOperator = read.elementType();

            if(currentOperator != null && !isValidOperatorTransition(currentOperator, readOperator))
                continue;

            if(operator == null)
            {
                operator = readOperator;
            }
            else if(operator != readOperator)
            {
                matched = false;
                break;
            }
        }

        if(matched)
            return operator;

        // follow rules to find consensus
        // if all reads are within the simplex tail at this location, and over 50% of them have the same cigar element, use that one

        // look at the reads that are duplex in this region. Classify the Cigar element within each duplex read as high qual or not:
        //  - For an element consuming read bases (e.g. M, I), it's high qual if the qual of the base itself is above the duplex mismatch qual = 1
        //  - For an element not consuming read bases (e.g. D), it’s always high qual (this is because a base is only deleted in the read
        //  if both R1 and R2 contain the deletion)

        // If over 50% of the reads that are duplex in this region contain the same Cigar element, with high quality, use that one
        // Otherwise, choose the most ref-like Cigar possible (e.g. strip I and D cigar elements)
        Map<String,Integer> consensusTypeCounts = Maps.newHashMap();
        int duplexCount = 0;
        int validReadCount = 0;
        boolean hasAligned = false;

        for(ReadParseState read : readStates)
        {
            if(read.beforeUnclippedPosition(refPosition))
                continue;

            CigarOperator readOperator = read.elementType();

            if(currentOperator != null && !isValidOperatorTransition(currentOperator, readOperator))
                continue;

            if(inAlignedBases && readOperator == S)
            {
                // ignore initial soft-clips once the aligned section has started
                if(isForwardConsensus && read.cigarIndex() == 0)
                    continue;
                else if(!isForwardConsensus && read.cigarIndex() > 0)
                    continue;
            }

            ++validReadCount;

            hasAligned |= readOperator == M;

            SbxQualType qualType;

            if(read.baseQual() == RAW_SIMPLEX_QUAL)
            {
                qualType = SbxQualType.SIMPLEX;
            }
            else
            {
                ++duplexCount;

                if(readOperator.consumesReadBases())
                    qualType = read.baseQual() > SBX_DUPLEX_MISMATCH_QUAL ? SbxQualType.DUPLEX_HIGH_QUAL : SbxQualType.DUPLEX_LOW_QUAL;
                else
                    qualType = SbxQualType.DUPLEX_HIGH_QUAL;
            }

            String key = formQualTypeOperatorKey(qualType, readOperator);
            Integer count = consensusTypeCounts.get(key);
            consensusTypeCounts.put(key, count != null ? count + 1 : 1);
        }

        CigarOperator consensusOperator = null;
        if(duplexCount == 0)
        {
            double minCount = validReadCount * 0.5;
            consensusOperator = getMajorityOperatorByType(consensusTypeCounts, SbxQualType.SIMPLEX, minCount);
        }
        else
        {
            double minCount = duplexCount * 0.5;
            consensusOperator = getMajorityOperatorByType(consensusTypeCounts, SbxQualType.DUPLEX_HIGH_QUAL, minCount);
        }

        if(consensusOperator != null)
            return consensusOperator;

        // favour an M over an indel
        if(hasAligned)
            return M;

        // take most common
        return getMajorityOperatorByType(consensusTypeCounts, null, 0);
    }

    private static CigarOperator getMajorityOperatorByType(
            final Map<String,Integer> consensusTypeCounts, @Nullable final SbxQualType requiredQualType, final double minCount)
    {
        CigarOperator maxOperator = null;
        int maxCount = 0;
        for(Map.Entry<String,Integer> entry : consensusTypeCounts.entrySet())
        {
            SbxQualType qualType = extractQualType(entry.getKey());

            if(requiredQualType == null || requiredQualType == qualType)
            {
                CigarOperator operator = extractCigarOperator(entry.getKey());
                int count = entry.getValue();

                if(maxCount == 0 || count > maxCount)
                {
                    maxCount = count;
                    maxOperator = operator;
                }
                else if(count == maxCount && useFirstCigarOperator(operator, maxOperator))
                {
                    maxCount = count;
                    maxOperator = operator;
                }
            }
        }

        return maxCount > minCount ? maxOperator : null;
    }

    private static boolean useFirstCigarOperator(final CigarOperator first, final CigarOperator second)
    {
        return first.ordinal() <= second.ordinal();
    }

    private static boolean isValidOperatorTransition(final CigarOperator current, final CigarOperator proposed)
    {
        if(current == proposed || current == M)
            return true;

        // all others must transition to aligned
        return proposed == M;
    }

    private enum SbxQualType
    {
        SIMPLEX,
        DUPLEX_LOW_QUAL,
        DUPLEX_HIGH_QUAL;
    }

    public static CigarOperator extractCigarOperator(final String key)
    {
        String[] components = key.split("_", 2);
        return CigarOperator.valueOf(components[0]);
    }

    public static SbxQualType extractQualType(final String key)
    {
        String[] components = key.split("_", 2);
        return SbxQualType.valueOf(components[1]);
    }

    public static String formQualTypeOperatorKey(final SbxQualType qualType, final CigarOperator operator)
    {
        return format("%s_%s", operator, qualType);
    }

}
