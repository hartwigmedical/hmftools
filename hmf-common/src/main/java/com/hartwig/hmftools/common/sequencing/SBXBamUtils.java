package com.hartwig.hmftools.common.sequencing;

import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftHardClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.replaceXwithM;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightHardClipLength;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;

import static htsjdk.samtools.CigarOperator.H;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;
import static htsjdk.samtools.SAMUtils.phredToFastq;

import java.util.Collections;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class SBXBamUtils
{
    public static byte INVALID_BASE_QUAL = -1;
    public static int INVALID_POSITION = -1;

    public static int DUPLEX_QUAL = 93;
    public static int SIMPLEX_QUAL = 18;
    public static int DUPLEX_ERROR_QUAL = 0;

    public static int BWA_MATCH_SCORE = 1;
    public static int BWA_MISMATCH_PENALTY = 4;
    public static int BWA_GAP_OPEN_PENALTY = 6;
    public static int BWA_GAP_EXTEND_PENALTY = 1;

    public static float ALIGNMENT_SCORED_DIFF_TO_MAPQ_DIFF = 1.0f;
    public static int MAX_MAPQ = 60;

    @Nullable
    private static String parseInt(final String s, int start)
    {
        if(start < 0 || start >= s.length())
            return null;

        if(s.charAt(start) < '0' || s.charAt(start) > '9')
            return null;

        StringBuilder intString = new StringBuilder();
        for(int i = start; i < s.length(); i++)
        {
            if(s.charAt(i) < '0' || s.charAt(i) > '9')
                break;

            intString.append(s.charAt(i));
        }

        return intString.toString();
    }

    public static List<Boolean> getDuplexIndels(final String ycTagStr)
    {
        List<Boolean> duplexIndels = Lists.newArrayList();
        String[] ycTagComponents = ycTagStr.split("-");

        int simplexHeadLength = Integer.parseInt(ycTagComponents[0]);
        String duplexRegion = ycTagComponents[1];
        for(int i = 0; i < simplexHeadLength; i++)
        {
            duplexIndels.add(false);
        }

        for(int i = 0; i < duplexRegion.length();)
        {
            String intString = parseInt(duplexRegion, i);
            if(intString != null)
            {
                int duplexMatchLength = Integer.parseInt(intString);
                for(int j = 0; j < duplexMatchLength; j++)
                {
                    duplexIndels.add(false);
                }
                i += intString.length();
                continue;
            }

            char code = duplexRegion.charAt(i);
            i++;
            switch(code)
            {
                case 'I':
                case 'L':
                case 'P':
                case 'Q':
                case 'J':
                case 'O':
                case 'X':
                case 'Z':
                    duplexIndels.add(true);
                    break;
                default:
                    duplexIndels.add(false);
            }
        }

        return duplexIndels;
    }

    public static void stripDuplexIndels(final SAMRecord record, final byte[] refBases, final int refStart)
    {
        if(record.getReadUnmappedFlag())
        {
            return;
        }

        String ycTagStr = record.getStringAttribute("YC");
        if(ycTagStr == null)
        {
            throw new IllegalArgumentException(format("Read must have YC tag: %s", record.getSAMString()));
        }

        boolean isForward = !record.getReadNegativeStrandFlag();

        List<Boolean> duplexIndels = getDuplexIndels(ycTagStr);
        if(!isForward)
        {
            Collections.reverse(duplexIndels);
        }

        int readLeftHardClipLength = leftHardClipLength(record);
        int readRightHardClipLength = rightHardClipLength(record);
        duplexIndels = duplexIndels.subList(readLeftHardClipLength, duplexIndels.size() - readRightHardClipLength);

        boolean hasDuplexIndels = false;
        for(boolean duplexIndel : duplexIndels)
        {
            hasDuplexIndels |= duplexIndel;
        }

        if(!hasDuplexIndels)
        {
            return;
        }

        List<AnnotatedBase> annotatedBases = getAnnotatedBases(record, duplexIndels);
        boolean readModified = processAnnotatedBases(annotatedBases, isForward, refBases, refStart);
        if(!readModified)
        {
            return;
        }

        int newAlignmentStart = INVALID_POSITION;
        StringBuilder newReadString = new StringBuilder();
        StringBuilder newBaseQualString = new StringBuilder();
        List<CigarOperator> newOps = Lists.newArrayList();
        for(AnnotatedBase annotatedBase : annotatedBases)
        {
            if(annotatedBase.deleted())
            {
                continue;
            }

            newOps.add(annotatedBase.Op);

            if(!annotatedBase.isReadBase())
            {
                continue;
            }

            if(newAlignmentStart == INVALID_POSITION && !annotatedBase.Op.isClipping())
            {
                newAlignmentStart = annotatedBase.RefPos;
            }

            newReadString.append((char) annotatedBase.ReadBase);
            newBaseQualString.append(phredToFastq(annotatedBase.qual()));
        }

        if(newAlignmentStart == INVALID_POSITION)
        {
            return;
        }

        List<CigarElement> newCigarElements = Lists.newArrayList();
        if(readLeftHardClipLength > 0)
        {
            newCigarElements.add(new CigarElement(readLeftHardClipLength, H));
        }

        CigarOperator currentOp = null;
        int currentLength = 0;
        for(CigarOperator op : newOps)
        {
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

            newCigarElements.add(new CigarElement(currentLength, currentOp));
            currentOp = op;
            currentLength = 1;
        }

        if(currentLength > 0)
        {
            newCigarElements.add(new CigarElement(currentLength, currentOp));
        }

        if(readRightHardClipLength > 0)
        {
            newCigarElements.add(new CigarElement(readRightHardClipLength, H));
        }

        int oldInsertGaps = 0;
        int oldTotalInsertLength = 0;
        for(CigarElement el : record.getCigar().getCigarElements())
        {
            if(el.getOperator() != I)
            {
                continue;
            }

            oldInsertGaps++;
            oldTotalInsertLength += el.getLength();
        }

        int newInsertGaps = 0;
        int newTotalInsertLength = 0;
        for(CigarElement el : newCigarElements)
        {
            if(el.getOperator() != I)
            {
                continue;
            }

            newInsertGaps++;
            newTotalInsertLength += el.getLength();
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

        int mapqDiff = round(ALIGNMENT_SCORED_DIFF_TO_MAPQ_DIFF * alignmentScoreDiff);
        if(mapqDiff != 0)
        {
            int currentMapq = record.getMappingQuality();
            int newMapq = currentMapq + mapqDiff;
            record.setMappingQuality(min(MAX_MAPQ, newMapq));
        }

        record.setReadString(newReadString.toString());
        record.setBaseQualityString(newBaseQualString.toString());
        record.setCigar(new Cigar(newCigarElements));
        record.setAlignmentStart(newAlignmentStart);
    }

    @VisibleForTesting
    public static class AnnotatedBase
    {
        public final int ReadIndex;
        public final int RefPos;
        public final CigarOperator Op;
        public final byte ReadBase;
        public final boolean IsDuplexIndel;

        private byte mQual;
        private boolean mDeleted;

        public AnnotatedBase(int readIndex, int refPos, final CigarOperator op, byte readBase, byte qual, boolean isDuplexIndel)
        {
            ReadIndex = readIndex;
            RefPos = refPos;
            Op = op;
            ReadBase = readBase;
            IsDuplexIndel = isDuplexIndel;

            mQual = qual;
            mDeleted = false;
        }

        public boolean isReadBase() { return Op.consumesReadBases(); }
        public boolean isRefBase() { return Op.consumesReferenceBases(); }

        public void deleteBase() { mDeleted = true; }
        public boolean deleted() { return mDeleted; }

        public byte qual() { return mQual; }
        public boolean setQual(byte qual) {
            if(mQual == qual)
                return false;

            mQual = qual;
            return true;
        }

        @Override
        public boolean equals(final Object o)
        {
            if(this == o)
            {
                return true;
            }

            if(!(o instanceof AnnotatedBase))
            {
                return false;
            }

            final AnnotatedBase that = (AnnotatedBase) o;
            return ReadIndex == that.ReadIndex && RefPos == that.RefPos && ReadBase == that.ReadBase && IsDuplexIndel == that.IsDuplexIndel
                    && mQual == that.mQual && mDeleted == that.mDeleted && Op == that.Op;
        }

        @Override
        public int hashCode()
        {
            int hash = ReadIndex;
            hash = 31 * hash + RefPos;
            hash = 31 * hash + Op.hashCode();
            hash = 31 * hash + (int) ReadBase;
            hash = 31 * hash + (IsDuplexIndel ? 1 : 0);
            hash = 31 * hash + (int) mQual;
            hash = 31 * hash + (mDeleted ? 1 : 0);

            return hash;
        }

        @Override
        public String toString()
        {
            return "AnnotatedBase{" +
                    "ReadIndex=" + ReadIndex +
                    ", RefPos=" + RefPos +
                    ", Op=" + Op +
                    ", ReadBase=" + ReadBase +
                    ", IsDuplexIndel=" + IsDuplexIndel +
                    ", mQual=" + mQual +
                    ", mDeleted=" + mDeleted +
                    '}';
        }
    }

    @VisibleForTesting
    public static List<AnnotatedBase> getAnnotatedBases(final SAMRecord record, final List<Boolean> duplexIndels)
    {
        byte[] quals = record.getBaseQualities();
        byte[] bases = record.getReadBases();
        int readIndex = 0;
        int refPos = record.getAlignmentStart() - leftSoftClipLength(record);
        List<AnnotatedBase> annotatedBases = Lists.newArrayList();
        for(CigarElement el : record.getCigar().getCigarElements())
        {
            if(el.getOperator() == H)
                continue;

            boolean isRead = el.getOperator().consumesReadBases();
            boolean isRef = el.getOperator() == S || el.getOperator().consumesReferenceBases();

            if(isRead && isRef)
            {
                for(int i = 0; i < el.getLength(); i++)
                {
                    annotatedBases.add(new AnnotatedBase(readIndex, refPos, el.getOperator(), bases[readIndex], quals[readIndex], duplexIndels.get(readIndex)));
                    readIndex++;
                    refPos++;
                }

                continue;
            }

            if(isRead)
            {
                for(int i = 0; i < el.getLength(); i++)
                {
                    annotatedBases.add(new AnnotatedBase(readIndex, refPos - 1, el.getOperator(), bases[readIndex], quals[readIndex], duplexIndels.get(readIndex)));
                    readIndex++;
                }

                continue;
            }

            if(isRef)
            {
                for(int i = 0; i < el.getLength(); i++)
                {
                    annotatedBases.add(new AnnotatedBase(readIndex - 1, refPos, el.getOperator(), INVALID_BASE_QUAL, INVALID_BASE_QUAL, false));
                    refPos++;
                }

                continue;
            }

            throw new IllegalStateException(format("CigarOperator(%s) in read(%s) with cigar(%s) consumes neither read or ref bases",
                    el.getOperator().name(), record.getReferenceName(), record.getCigarString()));
        }

        return annotatedBases;
    }

    @VisibleForTesting
    public static boolean processAnnotatedBases(final List<AnnotatedBase> annotatedBases, boolean isForward, final byte[] refBases,
            int refStart)
    {
        if(!isForward)
            Collections.reverse(annotatedBases);

        boolean readModified = false;
        for(int i = 0; i < annotatedBases.size();)
        {
            AnnotatedBase annotatedBase = annotatedBases.get(i);
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

                duplexIndelBases.append((char) annotatedBases.get(j).ReadBase);
                duplexIndelEnd = j;
            }

            // look backward to get start of the repeat
            int repeatStrIndex = duplexIndelBases.length() - 1;
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

            // set quals of first duplexIndelBases.length() of the repeat to zero
            int basesSeen = 0;
            for(int j = readRepeatStart; basesSeen < duplexIndelBases.length(); j++)
            {
                if(!annotatedBases.get(j).isReadBase())
                    continue;

                readModified |= annotatedBases.get(j).setQual((byte) 0);
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

            int refIndex = startRefRepeatPos - refStart;
            repeatStrIndex = duplexIndelBases.length() - 1;
            int refRepeatLength = 0;
            while(refIndex >= 0 && refIndex < refBases.length)
            {
                if(refBases[refIndex] != (byte) duplexIndelBases.charAt(repeatStrIndex))
                    break;

                refRepeatLength++;
                if(readRepeatLength <= refRepeatLength)
                    break;

                refIndex += isForward ? -1 : 1;

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
            int lastInsertIndex = INVALID_POSITION;
            for(int j = readRepeatStart; j <= duplexIndelEnd; j++)
            {
                if(annotatedBases.get(j).Op == I)
                    lastInsertIndex = j;
            }

            int trimCount = 0;
            if(lastInsertIndex != INVALID_POSITION)
            {
                for(int j = lastInsertIndex; j >= readRepeatStart && trimCount < trimLength; j--)
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

            // ensure last duplexIndelBases.length() untrimmed bases have qual 0
            basesSeen = 0;
            for(int j = duplexIndelEnd; j >= readRepeatStart && basesSeen < duplexIndelBases.length(); j--)
            {
                if(!annotatedBases.get(j).deleted() && annotatedBases.get(j).isReadBase())
                {
                    basesSeen++;
                    readModified |= annotatedBases.get(j).setQual((byte) 0);
                }
            }

            i = duplexIndelEnd + 1;
        }

        if(readModified)
        {
            if(!isForward)
                Collections.reverse(annotatedBases);
        }

        return readModified;
    }

    public static void fillQualZeroMismatchesWithRef(final SAMRecord record, final byte[] refBases, int refStart)
    {
        if(record.getReadUnmappedFlag())
            return;

        replaceXwithM(record);

        int refPos = record.getAlignmentStart();
        int readIndex = 0;
        byte[] readBases = record.getReadBases();
        byte[] baseQuals = record.getBaseQualities();
        int nmDiff = 0;
        int alignmentScoreDiff = 0;
        for(CigarElement el : record.getCigar().getCigarElements())
        {
            if(el.getOperator() == M)
            {
                for(int i = 0; i < el.getLength(); i++)
                {
                    int refBasesIndex = refPos + i - refStart;
                    if(refBasesIndex < 0)
                    {
                        i += -refBasesIndex - 1;
                        continue;
                    }

                    if(refBasesIndex >= refBases.length)
                        break;

                    int readBaseIndex = readIndex + i;
                    if(baseQuals[readBaseIndex] > 0)
                        continue;

                    baseQuals[readBaseIndex] = 1;

                    byte refBase = refBases[refBasesIndex];
                    byte readBase = readBases[readBaseIndex];
                    if(refBase != readBase)
                    {
                        readBases[readBaseIndex] = refBase;
                        nmDiff--;
                        alignmentScoreDiff += BWA_MISMATCH_PENALTY + BWA_MATCH_SCORE;
                    }
                }
            }

            if(el.getOperator().consumesReadBases())
                readIndex += el.getLength();

            if(el.getOperator().consumesReferenceBases())
                refPos += el.getLength();
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

        int mapqDiff = round(ALIGNMENT_SCORED_DIFF_TO_MAPQ_DIFF * alignmentScoreDiff);
        if(mapqDiff != 0)
        {
            int currentMapq = record.getMappingQuality();
            int newMapq = currentMapq + mapqDiff;
            record.setMappingQuality(min(MAX_MAPQ, newMapq));
        }
    }
}
