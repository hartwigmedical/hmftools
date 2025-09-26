package com.hartwig.hmftools.redux.consensus;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_GAP_EXTEND_PENALTY;
import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_GAP_OPEN_PENALTY;
import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_MATCH_SCORE;
import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_MISMATCH_PENALTY;
import static com.hartwig.hmftools.common.bam.CigarUtils.checkLeftAlignment;
import static com.hartwig.hmftools.common.bam.CigarUtils.leftHardClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightHardClipLength;
import static com.hartwig.hmftools.common.bam.ConsensusType.DUAL;
import static com.hartwig.hmftools.common.bam.ConsensusType.NONE;
import static com.hartwig.hmftools.common.bam.ConsensusType.SINGLE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_TYPE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.extractConsensusType;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.readToString;
import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_BYTES;
import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_N_BYTE;
import static com.hartwig.hmftools.common.codon.Nucleotides.baseIndex;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.RAW_DUPLEX_MISMATCH_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.RAW_DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_ADJACENT_1_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_ADJACENT_1_QUAL_REF_MATCH;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_ADJACENT_2_3_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_MISMATCH_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_READ_INDEX_TAG;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_SIMPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_YC_TAG;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.RAW_SIMPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.getDuplexIndelIndices;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ReduxConstants.SBX_CONSENSUS_BASE_THRESHOLD;
import static com.hartwig.hmftools.redux.consensus.BaseBuilder.INVALID_POSITION;
import static com.hartwig.hmftools.redux.consensus.BaseQualPair.NO_BASE;
import static com.hartwig.hmftools.redux.consensus.CommonUtils.findMostCommonBase;
import static com.hartwig.hmftools.redux.consensus.CommonUtils.findMostCommonBaseCount;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_FAIL;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_MATCH;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_MISMATCH;
import static com.hartwig.hmftools.redux.consensus.IlluminaRoutines.isDualStrandAndIsFirstInPair;
import static com.hartwig.hmftools.redux.consensus.IndelConsensusReads.deleteOrSplit;

import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;
import static htsjdk.samtools.CigarOperator.X;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.redux.ReduxConfig;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public final class SbxRoutines
{
    // internal SBX base qual values to maintain knowledge of simplex vs duplex mismatches
    protected static final byte DUPLEX_NO_CONSENSUS_QUAL = 3;
    protected static final byte SIMPLEX_NO_CONSENSUS_QUAL = 2;

    public static boolean SBX_HOMOPOLYMER_5_PRIME_LOW_BASE_QUAL = true; // currently unused but keep for now

    public static void stripDuplexIndels(final SAMRecord record)
    {
        if(record.getReadUnmappedFlag())
            return;

        // String chromosome = record.getReferenceName();
        String ycTagStr = record.getStringAttribute(SBX_YC_TAG);
        if(ycTagStr == null)
        {
            throw new IllegalArgumentException(format("read missing %s tag: %s", SBX_YC_TAG, readToString(record)));
        }

        List<Integer> duplexIndelIndices = getDuplexIndelIndices(ycTagStr);

        if(duplexIndelIndices == null)
            return;

        // not expecting to see hard-clips but remove any if present
        if(leftHardClipLength(record) > 0 || rightHardClipLength(record) > 0)
        {
            RD_LOGGER.error("hard-clipped reads not supported in SBX: {}", readToString(record));
            System.exit(1);
        }

        if(duplexIndelIndices.isEmpty())
            return;

        if(record.getReadNegativeStrandFlag())
            duplexIndelIndices = reverseDuplexIndelIndices(duplexIndelIndices, record.getReadBases().length);

        SbxDuplexIndelBuilder builder = new SbxDuplexIndelBuilder(record, duplexIndelIndices);

        List<SbxDuplexIndel> duplexIndels = builder.duplexIndels();

        if(duplexIndels.isEmpty())
            return;

        int lastDeletedIndex = -1;
        int netStrippedBases = 0;

        for(SbxDuplexIndel duplexIndel : duplexIndels)
        {
            if(lastDeletedIndex < 0)
            {
                netStrippedBases += duplexIndel.deletedBaseCount();
            }
            else
            {
                int adustedDeletedIndexStart = max(lastDeletedIndex + 1, duplexIndel.DeletedIndelIndexStart);
                netStrippedBases += max(duplexIndel.DeletedIndelIndexEnd - adustedDeletedIndexStart + 1, 0);
            }

            lastDeletedIndex = duplexIndel.DeletedIndelIndexEnd;

        }

        int oldBaseLength = record.getReadBases().length;
        int newBaseLength = record.getReadBases().length - netStrippedBases;

        byte[] readBases = new byte[newBaseLength];
        byte[] readQuals = new byte[newBaseLength];

        List<CigarElement> oldCigarElements = record.getCigar().getCigarElements();
        List<CigarElement> newCigarElements = Lists.newArrayListWithCapacity(oldCigarElements.size());

        int duplexIndelIndex = 0;
        SbxDuplexIndel duplexIndel = duplexIndels.get(duplexIndelIndex);
        int indelTrimmedCount = 0;

        int oldCigarIndex = 0;
        int oldCigarElementIndex = 0;
        CigarElement oldElement = oldCigarElements.get(0);

        CigarOperator curentCigarOp = null;
        int currentCigarLength = 0;

        int oldInsertCount = 0;
        int oldInsertedBases = 0;

        int newReadIndex = 0;

        for(int oldReadIndex = 0; oldReadIndex < oldBaseLength; ++oldCigarElementIndex)
        {
            if(oldCigarElementIndex >= oldElement.getLength())
            {
                // move to next old element, register insert info
                oldCigarElementIndex = 0;
                ++oldCigarIndex;
                oldElement = oldCigarElements.get(oldCigarIndex);

                if(oldElement.getOperator() == I)
                {
                    ++oldInsertCount;
                    oldInsertedBases += oldElement.getLength();
                }
            }

            boolean isStrippedIndelBase = false;
            boolean isLowQualBase = false;

            int effectReadIndex = newReadIndex + indelTrimmedCount; // factoring out trimmed bases
            boolean withinDuplexIndelBounds = duplexIndel != null && duplexIndel.withinBounds (effectReadIndex);

            if(withinDuplexIndelBounds && (oldElement.getOperator().consumesReadBases()))
            {
                if(duplexIndel.isLowQualBase(effectReadIndex))
                {
                    isLowQualBase = true;
                }
                else if(duplexIndel.withinDeleteBounds(effectReadIndex))
                {
                    // skip over stripped indel bases
                    isStrippedIndelBase = true;
                    ++indelTrimmedCount;
                }

                if(effectReadIndex >= duplexIndel.DuplexIndelIndexEnd)
                {
                    ++duplexIndelIndex;

                    if(duplexIndelIndex < duplexIndels.size())
                    {
                        duplexIndel = duplexIndels.get(duplexIndelIndex);

                        // skip over stripped indel bases if in next (adjacent) duplex indel
                        if(!isStrippedIndelBase && duplexIndel.withinDeleteBounds(effectReadIndex))
                        {
                            isStrippedIndelBase = true;
                            ++indelTrimmedCount;
                        }
                    }
                    else
                    {
                        duplexIndel = null;
                    }
                }
            }

            /*
            int effectReadIndex = newReadIndex + indelTrimmedCount; // factoring out trimmed bases
            boolean withinDuplexIndelBounds = duplexIndel != null && duplexIndel.withinBounds (effectReadIndex);

            if(withinDuplexIndelBounds && (oldElement.getOperator().consumesReadBases()))
            {
                if(duplexIndel.isLowQualBase(effectReadIndex))
                {
                    isLowQualBase = true;
                }
                else if(duplexIndel.withinDeleteBounds(effectReadIndex))
                {
                    // skip over stripped indel bases
                    isStrippedIndelBase = true;
                    ++indelTrimmedCount;
                }

                if(effectReadIndex >= duplexIndel.DuplexIndelIndexEnd)
                {
                    ++duplexIndelIndex;
                    // totalDeletedBases += duplexIndel.deletedBaseCount();
                    totalDeletedBases += indelTrimmedCount; // the duplexIndel deleted base count may be out-of-date due to overlaps

                    if(duplexIndelIndex < duplexIndels.size())
                    {
                        duplexIndel = duplexIndels.get(duplexIndelIndex);
                        indelTrimmedCount = 0;

                        // factor in deleted bases from prior stripped indels - using totalDeletedBases
                        if(totalDeletedBases > 0)
                            duplexIndel = new SbxDuplexIndel(duplexIndel, -totalDeletedBases);

                        // skip over stripped indel bases if in next (adjacent) duplex indel
                        if(!isStrippedIndelBase && duplexIndel.withinDeleteBounds(newReadIndex))
                            isStrippedIndelBase = true;
                    }
                    else
                    {
                        duplexIndel = null;
                    }
                }
            }
            */

            if(!isStrippedIndelBase)
            {
                if(curentCigarOp != oldElement.getOperator())
                {
                    if(currentCigarLength > 0)
                        newCigarElements.add(new CigarElement(currentCigarLength, curentCigarOp));

                    curentCigarOp = oldElement.getOperator();
                    currentCigarLength = 1;
                }
                else
                {
                    ++currentCigarLength;
                }
            }

            if(oldElement.getOperator().consumesReadBases())
            {
                if(!isStrippedIndelBase)
                {
                    if(newReadIndex >= readBases.length || oldReadIndex >= record.getReadBases().length)
                        break;

                    readBases[newReadIndex] = record.getReadBases()[oldReadIndex];

                    byte existingQual = record.getBaseQualities()[oldReadIndex];
                    byte qual;

                    if(isLowQualBase)
                        qual = SBX_DUPLEX_MISMATCH_QUAL;
                    else if(existingQual == RAW_DUPLEX_MISMATCH_QUAL && withinDuplexIndelBounds)
                        qual = RAW_DUPLEX_QUAL; // restore to duplex value
                    else
                        qual = existingQual;

                    readQuals[newReadIndex] = qual;
                    ++newReadIndex;
                }

                ++oldReadIndex;
            }
        }

        // add last
        newCigarElements.add(new CigarElement(currentCigarLength, curentCigarOp));

        int newInsertCount = 0;
        int newInsertedBases = 0;

        for(CigarElement element : newCigarElements)
        {
            if(element.getOperator() == I)
            {
                ++newInsertCount;
                newInsertedBases += element.getLength();
            }
        }

        int nmDiff = newInsertedBases - oldInsertedBases;
        Integer oldNumMutations = record.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE);
        if(oldNumMutations != null && nmDiff != 0)
        {
            int newNumMutations = oldNumMutations + nmDiff;
            record.setAttribute(NUM_MUTATONS_ATTRIBUTE, newNumMutations);
        }

        int oldInsertAlignmentScore = -BWA_GAP_OPEN_PENALTY * oldInsertCount - BWA_GAP_EXTEND_PENALTY * oldInsertedBases;
        int newInsertAlignmentScore = -BWA_GAP_OPEN_PENALTY * newInsertCount - BWA_GAP_EXTEND_PENALTY * newInsertedBases;
        int alignmentScoreDiff = newInsertAlignmentScore - oldInsertAlignmentScore;
        Integer oldAlignmentScore = record.getIntegerAttribute(ALIGNMENT_SCORE_ATTRIBUTE);
        if(oldAlignmentScore != null && alignmentScoreDiff != 0)
        {
            int newAlignmentScore = oldAlignmentScore + alignmentScoreDiff;
            record.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, newAlignmentScore);
        }

        record.setReadBases(readBases);
        record.setBaseQualities(readQuals);

        checkLeftAlignment(newCigarElements, readBases);

        record.setCigar(new Cigar(newCigarElements));

        if(ReduxConfig.RunChecks)
        {
            ReadValidReason validReason = ReadValidReason.isValidRead(record);

            if(validReason != ReadValidReason.OK)
            {
                RD_LOGGER.debug("invalid read({}) reason({}) details: {}", record.getReadName(), validReason, readToString(record));
            }
        }
    }

    private static List<Integer> reverseDuplexIndelIndices(List<Integer> duplexIndelIndices, int readLength)
    {
        List<Integer> newIndices = Lists.newArrayListWithCapacity(duplexIndelIndices.size());

        for(Integer index : duplexIndelIndices)
        {
            newIndices.add(0, readLength - index - 1);
        }

        return newIndices;
    }

    public static BaseQualPair determineBaseAndQual(
            final byte[] locationBases, final byte[] locationQuals, final String chromosome, int position, final RefGenome refGenome)
    {
        // count bases by qual type and apply rules
        Map<Byte,int[]> baseCountsByQual = Maps.newHashMap();

        int lowQualCount = 0;
        int simplexCount = 0;
        int duplexCount = 0;

        for(int i = 0; i < locationBases.length; ++i)
        {
            byte qual = locationQuals[i];

            if(qual == RAW_SIMPLEX_QUAL)
            {
                ++simplexCount;
            }
            else if(qual == RAW_DUPLEX_QUAL)
            {
                ++duplexCount;
            }
            else
            {
                ++lowQualCount;
                continue;
            }

            int[] qualBaseCounts = baseCountsByQual.get(qual);

            if(qualBaseCounts == null)
            {
                qualBaseCounts = new int[DNA_BASE_BYTES.length];
                baseCountsByQual.put(qual, qualBaseCounts);
            }

            int baseIndex = baseIndex(locationBases[i]);

            if(baseIndex >= 0 && baseIndex < DNA_BASE_BYTES.length)
            {
                ++qualBaseCounts[baseIndex];
            }
        }

        byte refBase = chromosome != null && position != INVALID_POSITION ? refGenome.getRefBase(chromosome, position) : DNA_N_BYTE;

        if(baseCountsByQual.isEmpty()) // all reads have invalid qual
            return new BaseQualPair(refBase, SBX_DUPLEX_MISMATCH_QUAL);

        if(simplexCount > 0 && duplexCount == 0)
        {
            int[] baseCounts = baseCountsByQual.get(RAW_SIMPLEX_QUAL);
            int maxBaseCount = findMostCommonBaseCount(baseCounts);
            byte maxBase = findMostCommonBase(baseCounts, refBase, maxBaseCount);

            int totalReadCount = simplexCount + lowQualCount;

            if(maxBaseCount > SBX_CONSENSUS_BASE_THRESHOLD * totalReadCount)
                return new BaseQualPair(maxBase, RAW_SIMPLEX_QUAL);
            else
                return new BaseQualPair(refBase, SIMPLEX_NO_CONSENSUS_QUAL);
        }

        if(duplexCount > 0)
        {
            int[] baseCounts = baseCountsByQual.get(RAW_DUPLEX_QUAL);
            int maxBaseCount = findMostCommonBaseCount(baseCounts);
            byte maxBase = findMostCommonBase(baseCounts, refBase, maxBaseCount);

            if(maxBaseCount > SBX_CONSENSUS_BASE_THRESHOLD * (duplexCount + lowQualCount))
                return new BaseQualPair(maxBase, RAW_DUPLEX_QUAL);
            else
                return new BaseQualPair(refBase, DUPLEX_NO_CONSENSUS_QUAL);
        }

        return new BaseQualPair(refBase, SBX_DUPLEX_MISMATCH_QUAL);
    }

    public static int findMaxDuplexBaseIndex(final List<SAMRecord> reads)
    {
        int firstDuplexBaseIndex = -1;

        // simplex bases are on the 5' end of the read, so look from that end for the first duplex base
        boolean fromStart = !reads.get(0).getReadNegativeStrandFlag();

        int maxReadIndex = reads.stream().mapToInt(x -> x.getReadBases().length).max().orElse(0);
        int index = fromStart ? 0 : maxReadIndex - 1;

        while(index >= 0 && index < maxReadIndex)
        {
            for(SAMRecord read : reads)
            {
                if(index >= read.getReadBases().length)
                    continue;

                if(read.getBaseQualities()[index] != RAW_SIMPLEX_QUAL)
                {
                    return index;
                }
            }

            index += fromStart ? 1 : -1;
        }

        return firstDuplexBaseIndex;
    }

    public static void determineIndelConsensus(
            final IndelConsensusReads indelConsensusReads, final ConsensusState consensusState, final List<SAMRecord> reads)
    {
        List<ReadParseState> readStates = reads.stream().map(x -> new ReadParseState(x, consensusState.IsForward)).collect(Collectors.toList());

        int[] consensusStarts = determineConsensusStartPositions(readStates);
        int consensusUnclippedStart = consensusStarts[0];
        int consensusAlignmentStart = consensusStarts[1];

        List<CigarElement> cigarElements = determineConsensusCigar(
                readStates, consensusState.IsForward, consensusUnclippedStart, consensusAlignmentStart);

        int baseLength = cigarElements.stream().filter(x -> x.getOperator().consumesReadBases()).mapToInt(x -> x.getLength()).sum();
        consensusState.setBaseLength(baseLength);

        int baseIndex = consensusState.IsForward ? 0 : baseLength - 1;

        int cigarCount = cigarElements.size();
        int cigarIndex = consensusState.IsForward ? 0 : cigarCount - 1;
        boolean[] isFirstInPair = null; // unused

        readStates.forEach(x -> x.reset());

        while(cigarIndex >= 0 && cigarIndex < cigarCount)
        {
            CigarElement element = cigarElements.get(cigarIndex);

            // simplest scenario is where all reads agree about this next element
            indelConsensusReads.addElementBases(consensusState, readStates, element, baseIndex, false, isFirstInPair);

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

        consensusState.finaliseCigar();

        int consensusAlignmentEnd = consensusAlignmentStart
                + cigarElements.stream().filter(x -> x.getOperator().consumesReferenceBases()).mapToInt(x -> x.getLength()).sum() - 1;

        int consensusUnclippedEnd = consensusAlignmentEnd;

        if(cigarElements.get(cigarElements.size() - 1).getOperator() == S)
            consensusUnclippedEnd += cigarElements.get(cigarElements.size() - 1).getLength();

        consensusState.setBoundaries(consensusUnclippedStart, consensusUnclippedEnd, consensusAlignmentStart, consensusAlignmentEnd);

        if(consensusState.outcome() != INDEL_FAIL)
            consensusState.setOutcome(INDEL_MISMATCH);
    }

    private static int[] determineConsensusStartPositions(final List<ReadParseState> readStates)
    {
        int earliestUnclippedStart = -1;
        int earliestAlignedStart = -1;

        for(ReadParseState read : readStates)
        {
            int alignmentStart = read.Read.getAlignmentStart();
            if(earliestAlignedStart < 0 || alignmentStart < earliestAlignedStart)
                earliestAlignedStart = alignmentStart;

            int unclippedStart = alignmentStart - leftSoftClipLength(read.Read.getCigar().getCigarElements());

            if(earliestUnclippedStart < 0 || unclippedStart < earliestUnclippedStart)
                earliestUnclippedStart = unclippedStart;
        }

        return new int[] { earliestUnclippedStart, earliestAlignedStart };
    }

    private static List<CigarElement> determineConsensusCigar(
            final List<ReadParseState> readStates, final boolean isForwardConsensus,
            int consensusUnclippedStart, int consensusAlignmentStart)
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

        /* need to handle duplicate jitter and then decide how to compare S vs M
        if(consensusUnclippedStart < consensusAlignmentStart)
        {
            cigarElements.add(new CigarElement(consensusAlignmentStart - consensusUnclippedStart, S));
        }
        */

        List<ReadParseState> activeReadStates = Lists.newArrayList(readStates);

        int elementLength = 0;
        CigarOperator elementOperator = null;

        while(true)
        {
            CigarOperator nextOperator = determineConsensusCigarOperator(activeReadStates);

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
            }

            ++elementLength;

            // move reads ahead depending on whether they agree with the consensus or not
            for(ReadParseState read : activeReadStates)
            {
                // check for element type differences:

                // first skip past any insert if the selected element is aligned
                if(nextOperator == M && read.elementType() == I)
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
                    if(nextOperator == M)
                    {
                        if(deleteOrSplit(read.elementType()) || read.elementType() == S)
                        {
                            //
                        }
                        else if(read.elementType() == I)
                        {
                            return Collections.emptyList();
                        }
                    }
                    else if(nextOperator == I)
                    {
                        if(read.elementType() == M || deleteOrSplit(read.elementType()))
                        {
                            moveNext = false;
                        }
                    }
                }

                if(moveNext)
                    read.moveNextBase();
            }

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

        // add the last
        if(isForwardConsensus)
            cigarElements.add(new CigarElement(elementLength, elementOperator));
        else
            cigarElements.add(0, new CigarElement(elementLength, elementOperator));

        return cigarElements;
    }

    private static CigarOperator determineConsensusCigarOperator(final List<ReadParseState> readStates)
    {
        CigarOperator operator = null;
        boolean matched = true;

        for(ReadParseState read : readStates)
        {
            if(operator == null)
            {
                operator = read.elementType();
            }
            else if(operator != read.elementType())
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
        boolean hasAligned = false;

        for(ReadParseState read : readStates)
        {
            SbxQualType qualType;
            CigarOperator readOperator = read.elementType();
            hasAligned |= readOperator == M;

            if(read.currentBaseQual() == RAW_SIMPLEX_QUAL)
            {
                qualType = SbxQualType.SIMPLEX;
            }
            else
            {
                ++duplexCount;

                if(readOperator.consumesReadBases())
                    qualType = read.currentBaseQual() > SBX_DUPLEX_MISMATCH_QUAL ? SbxQualType.DUPLEX_HIGH_QUAL : SbxQualType.DUPLEX_LOW_QUAL;
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
            double minCount = readStates.size() * 0.5;
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
            }
        }

        return maxCount > minCount ? maxOperator : null;
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

    private static final int[] ADJACENT_INDICES = {-3, -2, -1, 1, 2, 3};

    public static void finaliseRead(final RefGenomeInterface refGenome, final SAMRecord record)
    {
        // first pass map to final base quals and note duplex mismatches requiring adjacent base adjustments
        List<Integer> duplexMismatchIndices = null;
        int readLength = record.getBaseQualities().length;
        int lastReadIndex = readLength - 1;
        String chromosome = record.getReferenceName();
        byte[] newBaseQuals = record.getBaseQualities();

        ConsensusType consensusType = extractConsensusType(record);
        Integer firstDuplexBaseIndex = null;
        int duplexRegionStart = -1;
        int duplexRegionEnd = -1;

        if(consensusType == NONE)
        {
            firstDuplexBaseIndex = findMaxDuplexBaseIndex(List.of(record));
        }
        else
        {
            firstDuplexBaseIndex = record.getIntegerAttribute(SBX_DUPLEX_READ_INDEX_TAG);
        }

        if(firstDuplexBaseIndex != null && firstDuplexBaseIndex >= 0)
        {
            // mark the whole read as DUAL if it has any duplex bases - downstream processes will then use the duplex base index to
            // determine if a base is single or dual for consensus classification
            if(consensusType == SINGLE)
                record.setAttribute(CONSENSUS_TYPE_ATTRIBUTE, DUAL.toString());

            if(record.getReadNegativeStrandFlag())
            {
                duplexRegionStart = 0;
                duplexRegionEnd = firstDuplexBaseIndex;
            }
            else
            {
                duplexRegionStart = firstDuplexBaseIndex;
                duplexRegionEnd = lastReadIndex;
            }
        }

        for(int i = 0; i < readLength; ++i)
        {
            // values to handle and convert:
            // 93 - convert to 40
            // 18 - convert to 27
            // 1 - leave as-is but adjust adjacent bases
            // 3 (DUPLEX_NO_CONSENSUS_QUAL) - convert 1 and adjust adjacent bases
            // 2 (SIMPLEX_NO_CONSENSUS_QUAL) - convert to 1

            byte qual = newBaseQuals[i];

            if(qual == DUPLEX_NO_CONSENSUS_QUAL || qual <= SBX_DUPLEX_MISMATCH_QUAL)
            {
                if(duplexMismatchIndices == null)
                    duplexMismatchIndices = Lists.newArrayList();

                duplexMismatchIndices.add(i);
            }

            switch(qual)
            {
                case RAW_DUPLEX_QUAL:
                    newBaseQuals[i] = SBX_DUPLEX_QUAL;

                    // fix instances where consensus has set simplex in a duplex region
                    if(duplexRegionStart >= 0 && (i < duplexRegionStart || i > duplexRegionEnd))
                        newBaseQuals[i] = SBX_DUPLEX_MISMATCH_QUAL;
                    break;

                case RAW_SIMPLEX_QUAL:
                    newBaseQuals[i] = SBX_SIMPLEX_QUAL;

                    // as above
                    if(duplexRegionStart >= 0 && i >= duplexRegionStart && i <= duplexRegionEnd)
                        newBaseQuals[i] = SBX_DUPLEX_MISMATCH_QUAL;

                    break;

                case SIMPLEX_NO_CONSENSUS_QUAL:
                    newBaseQuals[i] = SBX_DUPLEX_MISMATCH_QUAL;
                    break;

                case DUPLEX_NO_CONSENSUS_QUAL:
                case RAW_DUPLEX_MISMATCH_QUAL:
                default:
                    newBaseQuals[i] = SBX_DUPLEX_MISMATCH_QUAL;
                    break;
            }
        }

        List<CigarElement> cigarElements = record.getCigar().getCigarElements();
        boolean requiresCigarUpdate = false;

        Map<Integer,Byte> duplexMismatchRefBase = null; // ref base at all mis-match indices/locations

        if(!record.getReadUnmappedFlag())
        {
            int refPos = record.getAlignmentStart();
            int readIndex = 0;
            byte[] readBases = record.getReadBases();
            int nmDiff = 0;
            int alignmentScoreDiff = 0;

            for(int i = 0; i < cigarElements.size(); ++i)
            {
                CigarElement element = cigarElements.get(i);

                if(element.getOperator() == X)
                {
                    // replace with M
                    cigarElements = Lists.newArrayList(cigarElements);
                    requiresCigarUpdate = true;
                    element = new CigarElement(element.getLength(), M);
                    cigarElements.set(i, element);
                }

                if(element.getOperator() != M)
                {
                    if(element.getOperator().consumesReadBases())
                        readIndex += element.getLength();

                    if(element.getOperator().consumesReferenceBases())
                        refPos += element.getLength();

                    continue;
                }

                for(int j = 0; j < element.getLength(); j++, readIndex++, refPos++)
                {
                    if(readIndex >= newBaseQuals.length)
                        break;

                    if(newBaseQuals[readIndex] > SBX_DUPLEX_MISMATCH_QUAL)
                        continue;

                    if(refPos < 1)
                        continue;

                    if(duplexMismatchRefBase == null)
                        duplexMismatchRefBase = Maps.newHashMap();

                    byte refBase = refGenome.getBase(chromosome, refPos);

                    duplexMismatchRefBase.put(readIndex, refBase);

                    byte readBase = readBases[readIndex];
                    if(refBase == readBase)
                        continue;

                    readBases[readIndex] = refBase;
                    nmDiff--;
                    alignmentScoreDiff += BWA_MISMATCH_PENALTY + BWA_MATCH_SCORE;
                }
            }

            Integer oldNumMutations = record.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE);
            if(oldNumMutations != null && nmDiff != 0)
            {
                int newNumMutations = oldNumMutations + nmDiff;
                record.setAttribute(NUM_MUTATONS_ATTRIBUTE, newNumMutations);
            }

            Integer oldAlignmentScore = record.getIntegerAttribute(ALIGNMENT_SCORE_ATTRIBUTE);
            if(oldAlignmentScore != null && alignmentScoreDiff != 0)
            {
                int newAlignmentScore = oldAlignmentScore + alignmentScoreDiff;
                record.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, newAlignmentScore);
            }

            if(requiresCigarUpdate)
            {
                int index = 0;
                while(index < cigarElements.size() - 1)
                {
                    CigarElement element = cigarElements.get(index);

                    int nextIndex = index + 1;
                    while(nextIndex < cigarElements.size())
                    {
                        CigarElement nextElement = cigarElements.get(nextIndex);

                        if(element.getOperator() == M && nextElement.getOperator() == M)
                        {
                            element = new CigarElement(element.getLength() + nextElement.getLength(), M);
                            cigarElements.set(index, element);
                            cigarElements.remove(nextElement);
                        }
                        else
                        {
                            break;
                        }
                    }

                    ++index;
                }

                record.setCigar(new Cigar(cigarElements));
            }
        }

        if(duplexMismatchIndices != null)
        {
            // apply duplex adjacent error logic within the duplex region
            for(Integer readIndex : duplexMismatchIndices)
            {
                Byte misMatchRefBase = duplexMismatchRefBase != null ? duplexMismatchRefBase.get(readIndex) : null;

                for(int j = 0; j < ADJACENT_INDICES.length; ++j)
                {
                    byte adjustQual = abs(ADJACENT_INDICES[j]) == 1 ? SBX_DUPLEX_ADJACENT_1_QUAL : SBX_DUPLEX_ADJACENT_2_3_QUAL;
                    int adjustIndex = readIndex + ADJACENT_INDICES[j];

                    if(adjustIndex < 0 || adjustIndex > lastReadIndex)
                        continue;

                    // for the immediately adjacent base, if it matches the ref at the mismatch base, use a different adjusted qual
                    if(adjustQual == SBX_DUPLEX_ADJACENT_1_QUAL && misMatchRefBase != null)
                    {
                        if(record.getReadBases()[adjustIndex] == misMatchRefBase)
                            adjustQual = SBX_DUPLEX_ADJACENT_1_QUAL_REF_MATCH;
                    }

                    if(duplexRegionStart >= 0 && (adjustIndex < duplexRegionStart || adjustIndex > duplexRegionEnd))
                        continue;

                    byte indexQual = newBaseQuals[adjustIndex];

                    if(indexQual > adjustQual && indexQual != SBX_SIMPLEX_QUAL)
                        newBaseQuals[adjustIndex] = adjustQual;
                }
            }
        }

        record.setBaseQualities(newBaseQuals);
    }

    // unused methods for now
    public static int determineBaseMatchQual(int existingQual, byte newQual)
    {
        if(existingQual < 0)
            return newQual;

        if(existingQual == newQual)
            return existingQual;

        if(newQual == RAW_DUPLEX_QUAL || existingQual == RAW_DUPLEX_QUAL)
            return RAW_DUPLEX_QUAL;

        // ensure a mismatch duplex qual takes precedence over simple quals
        if(newQual <= SBX_DUPLEX_MISMATCH_QUAL || existingQual <= SBX_DUPLEX_MISMATCH_QUAL)
            return SBX_DUPLEX_MISMATCH_QUAL;

        return max(existingQual, newQual);
    }

    private static boolean isDuplexQual(final byte qual)
    {
        return qual == RAW_DUPLEX_QUAL || qual == DUPLEX_NO_CONSENSUS_QUAL || qual == SBX_DUPLEX_MISMATCH_QUAL;
    }

    private static boolean isSimplexQual(final byte qual)
    {
        return qual == RAW_SIMPLEX_QUAL || qual == SIMPLEX_NO_CONSENSUS_QUAL;
    }
}
