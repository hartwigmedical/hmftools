package com.hartwig.hmftools.redux.consensus;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_GAP_EXTEND_PENALTY;
import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_GAP_OPEN_PENALTY;
import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_MATCH_SCORE;
import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_MISMATCH_PENALTY;
import static com.hartwig.hmftools.common.bam.CigarUtils.checkLeftAlignment;
import static com.hartwig.hmftools.common.bam.CigarUtils.collapseCigarOps;
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
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.reverseDuplexIndelIndices;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ReduxConstants.SBX_CONSENSUS_BASE_THRESHOLD;
import static com.hartwig.hmftools.redux.consensus.BaseBuilder.INVALID_POSITION;
import static com.hartwig.hmftools.redux.consensus.SbxAnnotatedBase.INVALID_BASE;
import static com.hartwig.hmftools.redux.consensus.SbxDuplexIndelBuilder.checkSupplementaryCigar;

import static htsjdk.samtools.CigarOperator.H;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;
import static htsjdk.samtools.CigarOperator.X;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public final class SbxRoutines
{
    // internal SBX base qual values to maintain knowledge of simplex vs duplex mismatches
    protected static final byte DUPLEX_NO_CONSENSUS_QUAL = 3;
    protected static final byte SIMPLEX_NO_CONSENSUS_QUAL = 2;

    public static boolean SBX_HOMOPOLYMER_5_PRIME_LOW_BASE_QUAL = true;

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

        SbxDuplexIndelBuilder builder = new SbxDuplexIndelBuilder(record, duplexIndelIndices);

        List<SbxDuplexIndel> duplexIndels = builder.duplexIndels();

        if(duplexIndels.isEmpty())
            return;

        int strippedBases = duplexIndels.stream().mapToInt(x -> x.deletedBaseCount()).sum();
        int oldBaseLength = record.getReadBases().length;
        int newBaseLength = record.getReadBases().length - strippedBases;

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

        int totalDeletedBases = 0;

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

            if(withinDuplexIndelBounds && (oldElement.getOperator() == I || oldElement.getOperator() == M))
            {
                if(duplexIndel.isLowQualBase(effectReadIndex))
                {
                    isLowQualBase = true;
                }
                else if(effectReadIndex >= duplexIndel.DeletedIndelIndexStart && effectReadIndex <= duplexIndel.DeletedIndelIndexEnd)
                {
                    // skip over stripped indel bases
                    isStrippedIndelBase = true;
                    ++indelTrimmedCount;
                }

                if(effectReadIndex >= duplexIndel.DuplexIndelIndexEnd)
                {
                    ++duplexIndelIndex;
                    totalDeletedBases += duplexIndel.deletedBaseCount();

                    if(duplexIndelIndex < duplexIndels.size())
                    {
                        duplexIndel = duplexIndels.get(duplexIndelIndex);
                        indelTrimmedCount = 0;

                        // factor in deleted bases from prior stripped indels - using totalDeletedBases
                        if(totalDeletedBases > 0)
                            duplexIndel = new SbxDuplexIndel(duplexIndel, totalDeletedBases);
                    }
                    else
                    {
                        duplexIndel = null;
                    }
                }
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
    }

    public static void stripDuplexIndelsOld(final RefGenomeInterface refGenome, final SAMRecord record)
    {
        if(record.getReadUnmappedFlag())
            return;

        String chromosome = record.getReferenceName();
        String ycTagStr = record.getStringAttribute(SBX_YC_TAG);
        if(ycTagStr == null)
        {
            throw new IllegalArgumentException(format("read missing %s tag: %s", SBX_YC_TAG, readToString(record)));
        }

        boolean isForward = !record.getReadNegativeStrandFlag();

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

        if(!isForward)
            reverseDuplexIndelIndices(duplexIndelIndices, record.getReadBases().length);

        List<SbxAnnotatedBase> annotatedBases = getAnnotatedBases(record, duplexIndelIndices);
        boolean readModified = processAnnotatedBases(refGenome, chromosome, annotatedBases, isForward);

        if(!readModified)
            return;

        int newAlignmentStart = INVALID_POSITION;
        List<CigarOperator> newOps = Lists.newArrayList();

        byte[] readBases = record.getReadBases();
        byte[] readQuals = record.getBaseQualities();
        boolean replaceReadBaseQuals = false;

        int newBaseLength = (int)annotatedBases.stream().filter(x -> !x.deleted() && x.isReadBase()).count();

        if(newBaseLength != record.getReadBases().length)
        {
            readBases = new byte[newBaseLength];
            readQuals = new byte[newBaseLength];
            replaceReadBaseQuals = true;
        }

        int readIndex = 0;
        for(SbxAnnotatedBase annotatedBase : annotatedBases)
        {
            if(annotatedBase.deleted())
                continue;

            CigarOperator cigarOperator = annotatedBase.originalOperator();
            newOps.add(cigarOperator);

            if(!cigarOperator.consumesReadBases())
                continue;

            if(newAlignmentStart == INVALID_POSITION && !cigarOperator.isClipping())
            {
                newAlignmentStart = annotatedBase.RefPos;
            }

            if(readIndex >= readBases.length)
                break;

            readBases[readIndex] = annotatedBase.ReadBase;
            readQuals[readIndex] = annotatedBase.qual();
            ++readIndex;
        }

        if(newAlignmentStart == INVALID_POSITION)
            return;

        List<CigarElement> newCigarElements = Lists.newArrayList();
        newCigarElements.addAll(collapseCigarOps(newOps));

        int oldInsertGaps = 0;
        int oldTotalInsertLength = 0;
        for(CigarElement element : record.getCigar().getCigarElements())
        {
            if(element.getOperator() != I)
                continue;

            oldInsertGaps++;
            oldTotalInsertLength += element.getLength();
        }

        int newInsertGaps = 0;
        int newTotalInsertLength = 0;
        for(CigarElement element : newCigarElements)
        {
            if(element.getOperator() != I)
                continue;

            newInsertGaps++;
            newTotalInsertLength += element.getLength();
        }

        int nmDiff = newTotalInsertLength - oldTotalInsertLength;
        Integer oldNumMutations = record.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE);
        if(oldNumMutations != null && nmDiff != 0)
        {
            int newNumMutations = oldNumMutations + nmDiff;
            record.setAttribute(NUM_MUTATONS_ATTRIBUTE, newNumMutations);
        }

        int oldInsertAlignmentScore = -BWA_GAP_OPEN_PENALTY * oldInsertGaps - BWA_GAP_EXTEND_PENALTY * oldTotalInsertLength;
        int newInsertAlignmentScore = -BWA_GAP_OPEN_PENALTY * newInsertGaps - BWA_GAP_EXTEND_PENALTY * newTotalInsertLength;
        int alignmentScoreDiff = newInsertAlignmentScore - oldInsertAlignmentScore;
        Integer oldAlignmentScore = record.getIntegerAttribute(ALIGNMENT_SCORE_ATTRIBUTE);
        if(oldAlignmentScore != null && alignmentScoreDiff != 0)
        {
            int newAlignmentScore = oldAlignmentScore + alignmentScoreDiff;
            record.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, newAlignmentScore);
        }

        if(replaceReadBaseQuals)
        {
            record.setReadBases(readBases);
            record.setBaseQualities(readQuals);
        }

        record.setCigar(new Cigar(newCigarElements));
        record.setAlignmentStart(newAlignmentStart);
    }

    @VisibleForTesting
    public static List<SbxAnnotatedBase> getAnnotatedBases(final SAMRecord record, final List<Integer> duplexIndelIndices)
    {
        List<SbxAnnotatedBase> annotatedBases = Lists.newArrayList();

        // use any supplementary data to set aligned bases to help with indel stripping
        SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(record);

        if(suppData != null && !HumanChromosome.contains(suppData.Chromosome))
            suppData = null;

        final List<CigarElement> readCigarElements = record.getCigar().getCigarElements();
        List<CigarElement> remoteSuppElements;
        List<CigarElement> localCigarElements;

        int leftSoftClipLength = readCigarElements.get(0).getOperator() == S ? readCigarElements.get(0).getLength() : 0;
        int lastElement = readCigarElements.size() - 1;
        int rightSoftClipLength = readCigarElements.get(lastElement).getOperator() == S ? readCigarElements.get(lastElement).getLength() : 0;
        boolean useLeftClipSuppData = false;

        if(suppData != null)
        {
            List<CigarElement> suppElements = CigarUtils.cigarElementsFromStr(suppData.Cigar);

            if(suppData.isForwardOrient() == record.getReadNegativeStrandFlag())
                Collections.reverse(suppElements);

            if(leftSoftClipLength > rightSoftClipLength)
            {
                remoteSuppElements = checkSupplementaryCigar(suppElements, leftSoftClipLength, true);
                localCigarElements = readCigarElements.subList(1, lastElement + 1); // remove soft-clip from original
                useLeftClipSuppData = true;
            }
            else
            {
                localCigarElements = readCigarElements.subList(0, lastElement);
                remoteSuppElements = checkSupplementaryCigar(suppElements, rightSoftClipLength, false);
            }
        }
        else
        {
            localCigarElements = readCigarElements;
            remoteSuppElements = Collections.emptyList();
        }

        byte[] quals = record.getBaseQualities();
        byte[] bases = record.getReadBases();
        int readIndex = 0;
        boolean useLocalAlignments = true;

        for(int e = 0; e <= 1; ++e)
        {
            List<CigarElement> cigarElements;

            if(e == 0)
            {
                if(useLeftClipSuppData)
                {
                    // process supplementary-data aligned bases first
                    useLocalAlignments = false;
                }
            }
            else
            {
                if(suppData == null)
                    break;

                useLocalAlignments = !useLocalAlignments;
            }

            int refPos;
            if(useLocalAlignments)
            {
                refPos = record.getAlignmentStart();

                if(!useLeftClipSuppData)
                    refPos -= leftSoftClipLength;

                cigarElements = localCigarElements;
            }
            else
            {
                refPos = suppData.Position;
                cigarElements = remoteSuppElements;
            }

            for(CigarElement element : cigarElements)
            {
                if(element.getOperator() == H)
                    continue;

                boolean isRead = element.getOperator().consumesReadBases();
                boolean isRef = element.getOperator() == S || element.getOperator().consumesReferenceBases();

                if(isRead && isRef) // soft-clip or M
                {
                    for(int i = 0; i < element.getLength(); i++)
                    {
                        annotatedBases.add(new SbxAnnotatedBase(
                                readIndex, refPos, element.getOperator(), bases[readIndex], quals[readIndex],
                                duplexIndelIndices.contains(readIndex)));
                        readIndex++;
                        refPos++;
                    }

                    continue;
                }

                if(isRead) // an insert
                {
                    for(int i = 0; i < element.getLength(); i++)
                    {
                        annotatedBases.add(new SbxAnnotatedBase(
                                readIndex, refPos - 1, element.getOperator(), bases[readIndex], quals[readIndex],
                                duplexIndelIndices.contains(readIndex)));
                        readIndex++;
                    }

                    continue;
                }

                if(isRef) // ie a delete
                {
                    for(int i = 0; i < element.getLength(); i++)
                    {
                        annotatedBases.add(new SbxAnnotatedBase(
                                readIndex - 1, refPos, element.getOperator(), INVALID_BASE, INVALID_BASE, false));
                        refPos++;
                    }
                }
            }
        }

        if(suppData != null)
        {
            int suppSoftIndexStart, suppSoftIndexEnd;

            if(leftSoftClipLength > rightSoftClipLength)
            {
                suppSoftIndexStart = 0;
                suppSoftIndexEnd = leftSoftClipLength - 1;
            }
            else
            {
                suppSoftIndexEnd = annotatedBases.size() - 1;
                suppSoftIndexStart = suppSoftIndexEnd - rightSoftClipLength + 1;
            }

            for(int i = suppSoftIndexStart; i <= suppSoftIndexEnd; ++i)
            {
                annotatedBases.get(i).setSoftClipSuppData(suppData.Chromosome);
            }
        }

        return annotatedBases;
    }

    @VisibleForTesting
    public static boolean processAnnotatedBases(
            final RefGenomeInterface refGenome, final String chromosome, final List<SbxAnnotatedBase> annotatedBases, boolean isForward)
    {
        int chromosomeLength = refGenome.getChromosomeLength(chromosome);
        Integer suppDataChromosomeLength = null;

        boolean requireReversedBases = SBX_HOMOPOLYMER_5_PRIME_LOW_BASE_QUAL ? isForward : !isForward;

        if(requireReversedBases)
            Collections.reverse(annotatedBases);

        boolean readModified = false;
        for(int i = 0; i < annotatedBases.size();)
        {
            SbxAnnotatedBase annotatedBase = annotatedBases.get(i);
            if(annotatedBase.Op == S || !annotatedBase.isReadBase())
            {
                i++;
                continue;
            }

            if(!annotatedBase.IsDuplexIndel)
            {
                i++;
                continue;
            }

            // look forward to get contiguous duplex indel bases
            StringBuilder duplexIndelBases = new StringBuilder();
            duplexIndelBases.append((char) annotatedBase.ReadBase);
            int duplexIndelEnd = i;
            for(int j = i + 1; j < annotatedBases.size(); j++)
            {
                if(!annotatedBases.get(j).isReadBase())
                    continue;

                if(!annotatedBases.get(j).IsDuplexIndel)
                    break;

                duplexIndelBases.append((char)annotatedBases.get(j).ReadBase);
                duplexIndelEnd = j;
            }

            // look backward to get start of the repeat
            int repeatStrIndex = duplexIndelBases.length() - 1;
            int duplexIndelMismatches = duplexIndelBases.length();
            int readRepeatLength = duplexIndelBases.length();
            int readRepeatStart = i;
            for(int j = i - 1; j >= 0; j--)
            {
                if(!annotatedBases.get(j).isReadBase())
                    continue;

                if(annotatedBases.get(j).ReadBase != (byte) duplexIndelBases.charAt(repeatStrIndex))
                    break;

                readRepeatStart = j;
                readRepeatLength++;

                repeatStrIndex--;
                if(repeatStrIndex < 0)
                    repeatStrIndex = duplexIndelBases.length() - 1;
            }

            // set quals of first duplexIndelBases.length() of the repeat to SBX_DUPLEX_MISMATCH_QUAL
            int basesSeen = 0;
            for(int j = readRepeatStart; basesSeen < duplexIndelBases.length(); j++)
            {
                if(!annotatedBases.get(j).isReadBase())
                    continue;

                readModified |= annotatedBases.get(j).setQual(SBX_DUPLEX_MISMATCH_QUAL);
                basesSeen++;
            }

            // find ref repeat length
            int startRefRepeatPos = INVALID_POSITION;
            for(int j = duplexIndelEnd; j >= 0; j--)
            {
                if(annotatedBases.get(j).isRefBase())
                {
                    startRefRepeatPos = annotatedBases.get(j).RefPos;
                    break;
                }
            }

            if(startRefRepeatPos == INVALID_POSITION)
            {
                i = duplexIndelEnd + 1;
                continue;
            }

            int refPos = startRefRepeatPos;
            repeatStrIndex = duplexIndelBases.length() - 1;
            int refRepeatLength = 0;

            String refPosChr = chromosome;
            int refChromosomeLength = chromosomeLength;

            if(annotatedBase.inSoftClipSuppData())
            {
                refPosChr = annotatedBase.softClipChromosome();

                if(!refPosChr.equals(chromosome))
                {
                    if(suppDataChromosomeLength == null)
                    {
                        suppDataChromosomeLength = refGenome.getChromosomeLength(refPosChr);
                    }

                    refChromosomeLength = suppDataChromosomeLength;
                }
            }

            while(refPos >= 1 && refPos <= refChromosomeLength)
            {
                byte refBase = refGenome.getBase(refPosChr, refPos);
                if(refBase != (byte) duplexIndelBases.charAt(repeatStrIndex))
                    break;

                refRepeatLength++;
                if(readRepeatLength <= refRepeatLength)
                    break;

                refPos += isForward ? -1 : 1;

                repeatStrIndex--;
                if(repeatStrIndex < 0)
                    repeatStrIndex = duplexIndelBases.length() - 1;
            }

            if(readRepeatLength <= refRepeatLength)
            {
                i = duplexIndelEnd + 1;
                continue;
            }

            // trim bases
            int trimLength = min(readRepeatLength - refRepeatLength, duplexIndelBases.length());
            int insertTrimStart = INVALID_POSITION;

            int jStart = requireReversedBases ? duplexIndelEnd : readRepeatStart;
            int jInc = requireReversedBases ? -1 : 1;

            for(int j = jStart; j >= readRepeatStart && j <= duplexIndelEnd; j += jInc)
            {
                if(annotatedBases.get(j).Op == I)
                {
                    insertTrimStart = j;
                    break;
                }
            }

            int trimCount = 0;
            if(insertTrimStart != INVALID_POSITION)
            {
                for(int j = insertTrimStart; j >= readRepeatStart && j <= duplexIndelEnd && trimCount < trimLength; j += jInc)
                {
                    if(!annotatedBases.get(j).isReadBase())
                        continue;

                    if(annotatedBases.get(j).Op != I)
                        break;

                    annotatedBases.get(j).deleteBase();
                    trimCount++;
                    readModified = true;
                }
            }

            // ensure first duplexIndelBases.length() left-aligned untrimmed bases have qual 0
            basesSeen = 0;
            for(int j = jStart; j >= readRepeatStart && j <= duplexIndelEnd && basesSeen < duplexIndelBases.length(); j += jInc)
            {
                if(!annotatedBases.get(j).deleted() && annotatedBases.get(j).isReadBase())
                {
                    basesSeen++;
                    readModified |= annotatedBases.get(j).setQual(SBX_DUPLEX_MISMATCH_QUAL);
                }
            }

            i = duplexIndelEnd + 1;
        }

        if(readModified && requireReversedBases)
            Collections.reverse(annotatedBases);

        return readModified;
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
