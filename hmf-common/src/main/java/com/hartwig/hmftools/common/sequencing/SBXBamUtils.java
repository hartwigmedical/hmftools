package com.hartwig.hmftools.common.sequencing;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftHardClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightHardClipLength;

import static htsjdk.samtools.CigarOperator.H;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.S;
import static htsjdk.samtools.SAMUtils.phredToFastq;

import java.util.Collections;
import java.util.List;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class SBXBamUtils
{
    public static int DUPLEX_QUAL = 93;
    public static int SIMPLEX_QUAL = 18;
    public static int DUPLEX_ERROR_QUAL = 0;

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

    public static List<Boolean> getDuplexIndels(final SAMRecord record)
    {
        List<Boolean> duplexIndels = Lists.newArrayList();

        String ycTagStr = record.getStringAttribute("YC");
        String[] ycTagComponents = ycTagStr.split("[-]");

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

    private static byte INVALID_BASE_QUAL = -1;
    private static int INVALID_POSITION = -1;

    private static class AnnotatedBase
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
    }

    public static void stripDuplexIndels(final SAMRecord record, final byte[] refBases, final int refStart)
    {
        if(record.getReadUnmappedFlag())
            return;

        boolean isForward = !record.getReadNegativeStrandFlag();

        List<Boolean> duplexIndels = getDuplexIndels(record);
        if(!isForward)
            Collections.reverse(duplexIndels);

        int readLeftHardClipLength = leftHardClipLength(record);
        int readRightHardClipLength = rightHardClipLength(record);
        duplexIndels = duplexIndels.subList(readLeftHardClipLength, duplexIndels.size() - readRightHardClipLength);

        boolean hasDuplexIndels = false;
        for(boolean duplexIndel : duplexIndels)
            hasDuplexIndels |= duplexIndel;

        if(!hasDuplexIndels)
            return;

        // annotate bases
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

        if(!readModified)
            return;

        if(!isForward)
            Collections.reverse(annotatedBases);

        int newAlignmentStart = INVALID_POSITION;
        StringBuilder newReadString = new StringBuilder();
        StringBuilder newBaseQualString = new StringBuilder();
        List<CigarOperator> newOps = Lists.newArrayList();
        for(AnnotatedBase annotatedBase : annotatedBases)
        {
            if(annotatedBase.deleted())
                continue;

            newOps.add(annotatedBase.Op);

            if(!annotatedBase.isReadBase())
                continue;

            if(newAlignmentStart == INVALID_POSITION && !annotatedBase.Op.isClipping())
                newAlignmentStart = annotatedBase.RefPos;

            newReadString.append((char) annotatedBase.ReadBase);
            newBaseQualString.append(phredToFastq(annotatedBase.qual()));
        }

        if(newAlignmentStart == INVALID_POSITION)
            return;

        List<CigarElement> newCigarElements = Lists.newArrayList();
        if(readLeftHardClipLength > 0)
            newCigarElements.add(new CigarElement(readLeftHardClipLength, H));

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
            newCigarElements.add(new CigarElement(currentLength, currentOp));

        if(readRightHardClipLength > 0)
            newCigarElements.add(new CigarElement(readRightHardClipLength, H));

        record.setReadString(newReadString.toString());
        record.setBaseQualityString(newBaseQualString.toString());
        record.setCigar(new Cigar(newCigarElements));
        record.setAlignmentStart(newAlignmentStart);
    }
}
