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
import static com.hartwig.hmftools.redux.consensus.ReadValidReason.checkIsValidRead;

import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;
import static htsjdk.samtools.CigarOperator.X;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.common.ReadInfo;

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

    public static void prepProcessRead(final SAMRecord record)
    {
        try
        {
            stripDuplexIndels(record);

            // ensure both sides of repeats are marked with duplex mismatch low-qual bases
            int readIndex = 0;
            byte[] readBases = record.getReadBases();
            byte[] baseQuals = record.getBaseQualities();

            for(int i = 0; i < record.getCigar().getCigarElements().size(); ++i)
            {
                CigarElement element = record.getCigar().getCigarElements().get(i);

                if(element.getOperator() == M)
                {
                    ensureMirroredMismatchRepeatQuals(readIndex, element.getLength(), readBases, baseQuals);
                }

                if(element.getOperator().consumesReadBases())
                    readIndex += element.getLength();
            }

            if(ReduxConfig.RunChecks)
                checkIsValidRead(record);
        }
        catch(Exception e)
        {
            RD_LOGGER.error("pre-process read({}) error: {}", ReadInfo.readToString(record), e.toString());
            e.printStackTrace();
        }
    }

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

        // convert inserts to soft-clips if required
        if(hasInvalidCigar(record))
        {
            List<CigarElement> fixedElements = Lists.newArrayList(record.getCigar().getCigarElements());
            correctInvalidCigar(fixedElements);
            record.setCigar(new Cigar(fixedElements));
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

        // gather unique deleted indices - these can be present in multiple duplex indels
        Set<Integer> netStrippedBases = Sets.newHashSet();

        for(SbxDuplexIndel duplexIndel : duplexIndels)
        {
            for(int i = duplexIndel.DeletedIndelIndexStart; i <= duplexIndel.DeletedIndelIndexEnd; ++i)
            {
                netStrippedBases.add(i);
            }
        }

        int oldBaseLength = record.getReadBases().length;
        int newBaseLength = record.getReadBases().length - netStrippedBases.size();

        byte[] readBases = new byte[newBaseLength];
        byte[] readQuals = new byte[newBaseLength];

        List<CigarElement> oldCigarElements = record.getCigar().getCigarElements();
        List<CigarElement> newCigarElements = Lists.newArrayListWithCapacity(oldCigarElements.size());

        int minDuplexIndelIndex = 0;
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
            boolean withinDuplexIndelBounds = false;

            int effectReadIndex = newReadIndex + indelTrimmedCount; // factoring out trimmed bases

            if(oldElement.getOperator().consumesReadBases())
            {
                int i = minDuplexIndelIndex;
                int firstActiveIndex = -1;

                for(; i < duplexIndels.size(); ++i)
                {
                    SbxDuplexIndel duplexIndel = duplexIndels.get(i);

                    if(effectReadIndex > duplexIndel.DuplexIndelIndexEnd)
                        continue;

                    if(firstActiveIndex < 0)
                        firstActiveIndex = i;

                    if(!duplexIndel.withinBounds(effectReadIndex))
                        continue;

                    withinDuplexIndelBounds = true;

                    if(duplexIndel.withinDeleteBounds(effectReadIndex))
                    {
                        // skip over stripped indel bases
                        isStrippedIndelBase = true;
                        break;
                    }
                    if(duplexIndel.isLowQualBase(effectReadIndex))
                    {
                        isLowQualBase = true;
                    }
                }

                if(firstActiveIndex > 0)
                    minDuplexIndelIndex = firstActiveIndex; // to avoid re-checking duplex indels which have been fully processed
            }

            if(isStrippedIndelBase)
            {
                ++indelTrimmedCount;
                isLowQualBase = false;
            }

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

        // convert any initial insert to soft-clip, as keeping with alignment expectations
        correctInvalidCigar(newCigarElements);

        checkLeftAlignment(newCigarElements, readBases);

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

    @VisibleForTesting
    private static boolean hasInvalidCigar(final SAMRecord record)
    {
        return hasInvalidCigar(record.getCigar().getCigarElements());
    }

    @VisibleForTesting
    public static boolean hasInvalidCigar(final List<CigarElement> cigarElements)
    {
        // looks for a specific SBX alignment issue where an insert comes before the first aligned section
        if(cigarElements.get(0).getOperator() == I)
            return true;

        if(cigarElements.size() > 1)
        {
            int lastIndex = cigarElements.size() - 1;

            if(cigarElements.get(lastIndex).getOperator() == I)
                return true;

            if(cigarElements.get(0).getOperator() == S && cigarElements.get(1).getOperator() == I)
                return true;

            if(cigarElements.get(lastIndex - 1).getOperator() == I && cigarElements.get(lastIndex).getOperator() == S)
                return true;
        }

        return false;
    }

    @VisibleForTesting
    public static void correctInvalidCigar(final List<CigarElement> cigarElements)
    {
        // no attempt is made at this stage to correct the alignment score or num of mutations
        int lastIndex = cigarElements.size() - 1;

        int i = 0;
        while(i < cigarElements.size())
        {
            CigarElement element = cigarElements.get(i);

            if(element.getOperator() == I)
            {
                if(i == 0 || i == lastIndex)
                {
                    cigarElements.set(i, new CigarElement(cigarElements.get(i).getLength(), S));
                }
                else if(i == 1 && cigarElements.get(0).getOperator() == S)
                {
                    // combine initial soft-clip and insert
                    cigarElements.remove(i);
                    --lastIndex;
                    cigarElements.set(0, new CigarElement(cigarElements.get(0).getLength() + element.getLength(), S));
                    continue;
                }
                else if(i == lastIndex - 1 && cigarElements.get(lastIndex).getOperator() == S)
                {
                    // combine initial soft-clip and insert
                    cigarElements.remove(i);
                    --lastIndex;
                    cigarElements.set(lastIndex, new CigarElement(cigarElements.get(lastIndex).getLength() + element.getLength(), S));
                    continue;
                }
            }

            ++i;
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

    @VisibleForTesting
    public static void ensureMirroredMismatchRepeatQuals(int readIndexStart, int elementLength, final byte[] bases, final byte[] baseQuals)
    {
        int readIndexEnd = readIndexStart + elementLength - 1;
        for(int readIndex = readIndexStart; readIndex <= readIndexEnd; readIndex++)
        {
            if(readIndex >= bases.length)
                break;

            if(baseQuals[readIndex] > SBX_DUPLEX_MISMATCH_QUAL)
                continue;

            int repeatLength = 1;
            boolean isHomopolymer = true;
            byte repeatBase = bases[readIndex];

            for(int i = readIndex + 1; i <= readIndexEnd; ++i)
            {
                if(baseQuals[i] > SBX_DUPLEX_MISMATCH_QUAL)
                    break;

                if(bases[i] != repeatBase)
                    isHomopolymer = false;

                ++repeatLength;
            }

            if(isHomopolymer)
                repeatLength = 1;

            byte[] repeatBases = new byte[repeatLength];

            for(int i = 0; i < repeatLength; ++i)
            {
                repeatBases[i] = bases[readIndex + i];
            }

            // search for the repeat as far as possible in both directions
            int mirroredIndexStart = -1;
            int maxRepeatShift = 0;

            for(int j = 0; j <= 1; ++j)
            {
                boolean searchUp = (j == 0);

                int repeatStart = readIndex + (searchUp ? repeatLength : -repeatLength);
                int lastRepeatStartIndex = readIndex;

                while(repeatStart >= 0 && repeatStart + repeatLength - 1 <= readIndexEnd)
                {
                    boolean matched = true;

                    for(int i = 0; i < repeatLength; ++i)
                    {
                        if(repeatBases[i] != bases[repeatStart + i])
                        {
                            matched = false;
                            break;
                        }
                    }

                    if(!matched)
                        break;

                    lastRepeatStartIndex = repeatStart;
                    repeatStart += searchUp ? repeatLength : -repeatLength;
                }

                if(lastRepeatStartIndex != readIndex)
                {
                    int repeatShift = abs(lastRepeatStartIndex - readIndex);

                    if(repeatShift > maxRepeatShift)
                    {
                        maxRepeatShift = repeatShift;
                        mirroredIndexStart = lastRepeatStartIndex;
                    }
                }
            }

            // mark the last repeat length bases as mismatch
            if(mirroredIndexStart >= 0)
            {
                for(int i = 0; i < repeatLength; ++i)
                {
                    baseQuals[mirroredIndexStart + i] = SBX_DUPLEX_MISMATCH_QUAL;
                }
            }

            if(mirroredIndexStart > readIndex)
                readIndex = mirroredIndexStart + repeatLength - 1; // move to end of the repeat if shifted right
        }
    }

    private static int findMostCommonBaseCount(final int[] baseCounts)
    {
        int maxCount = 0;

        for(int b = 0; b < DNA_BASE_BYTES.length; ++b)
        {
            maxCount = max(maxCount, baseCounts[b]);
        }

        return maxCount;
    }

    private static byte findMostCommonBase(final int[] baseCounts, final byte refBase, final int maxCount)
    {
        byte maxBase = DNA_N_BYTE;

        for(int b = 0; b < DNA_BASE_BYTES.length; ++b)
        {
            if(baseCounts[b] == maxCount)
            {
                if(DNA_BASE_BYTES[b] == refBase || maxBase == DNA_N_BYTE)
                    maxBase = DNA_BASE_BYTES[b];
            }
        }

        return maxBase;
    }
}
