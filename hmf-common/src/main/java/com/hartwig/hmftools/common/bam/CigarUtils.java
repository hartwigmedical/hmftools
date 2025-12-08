package com.hartwig.hmftools.common.bam;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.INVALID_READ_INDEX;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;

import static htsjdk.samtools.CigarOperator.H;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;
import static htsjdk.samtools.CigarOperator.X;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public final class CigarUtils
{
    public static List<CigarElement> cigarElementsFromStr(final String cigarStr)
    {
        if(cigarStr.equals(NO_CIGAR))
            return Collections.emptyList();

        List<CigarElement> cigarElements = Lists.newArrayList();

        int length = 0;
        for(int i = 0; i < cigarStr.length(); ++i)
        {
            char c = cigarStr.charAt(i);
            int digit = c - '0';

            if(digit >= 0 && digit <= 9)
            {
                length = length * 10 + digit;
            }
            else
            {
                CigarOperator operator = CigarOperator.characterToEnum(c);
                cigarElements.add(new CigarElement(length, operator));
                length = 0;
            }
        }

        return cigarElements;
    }

    public static Cigar cigarFromStr(final String cigarStr)
    {
        return new Cigar(cigarElementsFromStr(cigarStr));
    }

    public static String cigarElementsToStr(final List<CigarElement> elements)
    {
        return elements.stream().map(x -> format("%d%s", x.getLength(), x.getOperator())).collect(Collectors.joining());
    }

    public static int cigarBaseLength(final Cigar cigar)
    {
        return cigar.getCigarElements().stream().filter(x -> x.getOperator().consumesReadBases()).mapToInt(x -> x.getLength()).sum();
    }

    public static int cigarAlignedLength(final Cigar cigar)
    {
        return cigar.getCigarElements().stream().filter(x -> x.getOperator().consumesReferenceBases()).mapToInt(x -> x.getLength()).sum();
    }

    public static int calcCigarAlignedLength(final String cigarStr)
    {
        // a string-parsing version of the method above
        int baseLength = 0;
        int currentElementLength = 0;
        for(int i = 0; i < cigarStr.length(); ++i)
        {
            char c = cigarStr.charAt(i);
            boolean isAddItem = (c == 'D' || c == 'M' || c == 'N');

            if(isAddItem)
            {
                baseLength += currentElementLength;
                currentElementLength = 0;
                continue;
            }
            int digit = c - '0';
            if(digit >= 0 && digit <= 9)
            {
                currentElementLength = currentElementLength * 10 + digit;
            }
            else
            {
                // we are in one of the ignored items such as S, H, N etc
                currentElementLength = 0;
            }
        }

        return baseLength;
    }

    public static boolean leftSoftClipped(final SAMRecord record) { return leftSoftClipped(record.getCigar()); }
    public static boolean rightSoftClipped(final SAMRecord record) { return rightSoftClipped(record.getCigar()); }

    public static boolean leftSoftClipped(final Cigar cigar)
    {
        return cigar.getCigarElements().get(0).getOperator() == S;
    }

    public static boolean rightSoftClipped(final Cigar cigar)
    {
        return cigar.getCigarElements().get(cigar.getCigarElements().size() - 1).getOperator() == S;
    }

    public static int leftSoftClipLength(final SAMRecord record) { return leftSoftClipLength(record.getCigar()); }

    public static int rightSoftClipLength(final SAMRecord record) { return rightSoftClipLength(record.getCigar()); }

    public static int leftSoftClipLength(final Cigar cigar)
    {
        CigarElement firstElement = cigar.getFirstCigarElement();
        return (firstElement != null && firstElement.getOperator() == S) ? firstElement.getLength() : 0;
    }

    public static int leftSoftClipLength(final List<CigarElement> elements)
    {
        return elements.get(0).getOperator() == S ? elements.get(0).getLength() : 0;
    }

    public static int rightSoftClipLength(final Cigar cigar)
    {
        CigarElement lastElement = cigar.getLastCigarElement();
        return (lastElement != null && lastElement.getOperator() == S) ? lastElement.getLength() : 0;
    }

    public static int leftHardClipLength(final SAMRecord record) { return leftHardClipLength(record.getCigar()); }

    public static int rightHardClipLength(final SAMRecord record) { return rightHardClipLength(record.getCigar()); }

    public static int leftHardClipLength(final Cigar cigar)
    {
        CigarElement firstElement = cigar.getFirstCigarElement();
        return (firstElement != null && firstElement.getOperator() == H) ? firstElement.getLength() : 0;
    }

    public static int rightHardClipLength(final Cigar cigar)
    {
        CigarElement lastElement = cigar.getLastCigarElement();
        return (lastElement != null && lastElement.getOperator() == H) ? lastElement.getLength() : 0;
    }

    @Nullable
    public static String leftSoftClipBases(@NotNull final SAMRecord record)
    {
        int leftClip = leftSoftClipLength(record);
        if(leftClip == 0)
            return null;

        return record.getReadString().substring(0, leftClip);
    }

    @Nullable
    public static String rightSoftClipBases(@NotNull final SAMRecord record)
    {
        int rightClip = rightSoftClipLength(record);
        if(rightClip == 0)
            return null;

        return record.getReadString().substring(record.getReadString().length() - rightClip);
    }

    public static int getReadBoundaryPosition(
            final int readStart, final String cigarStr, final boolean getReadStart, boolean includeSoftClipped)
    {
        if(getReadStart && !includeSoftClipped) // avoid a parse for the basic case
            return readStart;

        // gets either the read start position or read end position, either with or without soft-clipped bases
        int currentPosition = readStart;
        int elementLength = 0;

        for(int i = 0; i < cigarStr.length(); ++i)
        {
            char c = cigarStr.charAt(i);
            boolean isAddItem = (c == 'D' || c == 'M' || c == 'S' || c == 'N');

            if(isAddItem)
            {
                if(getReadStart)
                {
                    // back out the left clip if present
                    return c == 'S' ? readStart - elementLength : readStart;
                }

                if(c == 'S' && readStart == currentPosition)
                {
                    // ignore left-clip since getting the read end position
                }
                else if(c == 'S' && !includeSoftClipped)
                {
                    break;
                }
                else
                {
                    currentPosition += elementLength;
                }

                elementLength = 0;
                continue;
            }

            int digit = c - '0';
            if(digit >= 0 && digit <= 9)
            {
                elementLength = elementLength * 10 + digit;
            }
            else
            {
                elementLength = 0;
            }
        }

        // always pointing to the start of the next element, so need to move back a base
        return currentPosition - 1;
    }

    public static int getReadIndexFromPosition(final int alignmentStart, final List<CigarElement> cigarElements, int position)
    {
        return getReadIndexFromPosition(alignmentStart, cigarElements, position, false, false);
    }

    public static int getReadIndexFromPosition(
            final int alignmentStart, final List<CigarElement> cigarElements, int position, boolean lastIfGapped, boolean inferFromSoftClip)
    {
        if(position < alignmentStart && !inferFromSoftClip)
            return INVALID_READ_INDEX;

        int refPosition = alignmentStart;
        int index = 0;

        for(CigarElement element : cigarElements)
        {
            if(element.getOperator() == S && refPosition == alignmentStart && inferFromSoftClip)
                refPosition -= element.getLength();

            boolean consumesRefBases = element.getOperator().consumesReferenceBases() || (element.getOperator() == S && inferFromSoftClip);

            if(consumesRefBases && refPosition + element.getLength() - 1 >= position)
            {
                if(element.getOperator().consumesReadBases())
                {
                    index += position - refPosition;
                }
                else
                {
                    if(lastIfGapped)
                        --index;
                    else
                        return INVALID_READ_INDEX;
                }

                return index;
            }

            if(consumesRefBases)
                refPosition += element.getLength();

            if(element.getOperator().consumesReadBases())
                index += element.getLength();
        }

        return INVALID_READ_INDEX;
    }

    public static int getPositionFromReadIndex(final int alignmentStart, final List<CigarElement> cigarElements, int readIndex)
    {
        return getPositionFromReadIndex(alignmentStart, cigarElements, readIndex, false, false)[0];
    }

    public static final int[] NO_POSITION_INFO = { NO_POSITION, 0 };

    public static int[] getPositionFromReadIndex(
            final int alignmentStart, final List<CigarElement> cigarElements, int readIndex, boolean lastIfInserted, boolean inferFromSoftClip)
    {
        if(readIndex < 0)
            return NO_POSITION_INFO;

        int refPosition = alignmentStart;
        int index = 0;

        for(CigarElement element : cigarElements)
        {
            if(element.getOperator() == S && refPosition == alignmentStart && inferFromSoftClip)
                refPosition -= element.getLength();

            boolean consumesRefBases = element.getOperator().consumesReferenceBases() || (element.getOperator() == S && inferFromSoftClip);

            if(element.getOperator().consumesReadBases() && index + element.getLength() - 1 >= readIndex)
            {
                int readIndexShift = 0;

                if(consumesRefBases)
                {
                    refPosition += readIndex - index;
                }
                else
                {
                    if(!lastIfInserted)
                        return NO_POSITION_INFO;

                    --refPosition;
                    readIndexShift = index - readIndex - 1;
                }

                return new int[] { refPosition, readIndexShift };
            }

            if(consumesRefBases)
                refPosition += element.getLength();

            if(element.getOperator().consumesReadBases())
                index += element.getLength();
        }

        return NO_POSITION_INFO;
    }

    public static List<CigarElement> collapseCigarOps(final List<CigarOperator> ops)
    {
        List<CigarElement> elems = Lists.newArrayList();
        CigarOperator currentOp = null;
        int currentLength = 0;
        for(CigarOperator op : ops)
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

            elems.add(new CigarElement(currentLength, currentOp));
            currentOp = op;
            currentLength = 1;
        }

        if(currentLength > 0)
            elems.add(new CigarElement(currentLength, currentOp));

        return elems;
    }

    public static boolean checkLeftAlignment(final List<CigarElement> cigarElements, final byte[] readBases)
    {
        // shift any right-aligned inserts back to the left
        boolean modified = false;
        boolean convertFirstInsert = false;

        int readIndex = 0;

        for(int i = 0; i < cigarElements.size() - 1; ++i)
        {
            CigarElement element = cigarElements.get(i);
            CigarElement nextElement = cigarElements.get(i + 1);
            CigarElement prevElement = i > 0 ? cigarElements.get(i - 1) : null;

            boolean checkLeftShift = prevElement != null && element.getOperator() == I
                    && prevElement.getOperator() == M && nextElement.getOperator() == M;

            if(checkLeftShift)
            {
                // cannot shift back through the previous element, but can convert it exactly to a soft-clip - see below
                if(element.getLength() > prevElement.getLength())
                    checkLeftShift = false;
                else if(element.getLength() == prevElement.getLength() && i != 1)
                    checkLeftShift = false;
            }

            if(checkLeftShift)
            {
                byte[] repeatBases = new byte[element.getLength()];
                boolean isHomopolymer = true;

                for(int j = 0; j < repeatBases.length; ++j)
                {
                    if(readIndex + j >= readBases.length) // invalid read bases if insert cannot extract corresponding bases
                        return false;

                    repeatBases[j] = readBases[readIndex + j];

                    if(j >= 1 && repeatBases[j] != repeatBases[0])
                        isHomopolymer = false;
                }

                if(isHomopolymer && repeatBases.length > 1)
                {
                    repeatBases = new byte[] { repeatBases[0] };
                }

                // search backwards through the previous aligned section
                int priorIndex = readIndex;
                int matchedRepeatCount = 0;

                while(true)
                {
                    priorIndex -= repeatBases.length;

                    if(priorIndex < 0)
                        break;

                    boolean matches = true;

                    for(int j = 0; j < repeatBases.length; ++j)
                    {
                        if(repeatBases[j] != readBases[priorIndex + j])
                        {
                            matches = false;
                            break;
                        }
                    }

                    if(!matches)
                        break;

                    ++matchedRepeatCount;

                    // cannot search back through the previous element
                    if((matchedRepeatCount + 1) * repeatBases.length > prevElement.getLength())
                        break;
                }

                if(matchedRepeatCount > 0)
                {
                    int alignmentShift = matchedRepeatCount * repeatBases.length;

                    if(alignmentShift == prevElement.getLength())
                    {
                        if(i == 1)
                        {
                            // the first element will be removed and the insert converted to a soft-clip
                            convertFirstInsert = true;
                        }
                        else
                        {
                            alignmentShift = (matchedRepeatCount - 1) * repeatBases.length;
                        }
                    }

                    prevElement = new CigarElement(max(prevElement.getLength() - alignmentShift, 1), prevElement.getOperator());
                    cigarElements.set(i - 1, prevElement);

                    nextElement = new CigarElement(nextElement.getLength() + alignmentShift, nextElement.getOperator());
                    cigarElements.set(i + 1, nextElement);

                    readIndex -= alignmentShift; // to move back to start of next aligned section

                    modified = true;
                }
            }

            if(element.getOperator().consumesReadBases())
                readIndex += element.getLength();
        }

        if(convertFirstInsert)
        {
            cigarElements.remove(0);
            CigarElement firstElement = cigarElements.get(0);
            cigarElements.set(0, new CigarElement(firstElement.getLength(), S));
        }

        return modified;
    }

    public static boolean hasValidCigar(final List<CigarElement> cigarElements)
    {
        if(cigarElements.isEmpty())
            return false;
        
        int cigarCount = cigarElements.size();

        if(cigarCount > 1)
        {
            CigarOperator lastOperator = cigarElements.get(0).getOperator();

            if(lastOperator != S && lastOperator != M)
                return false;

            for(int i = 1; i < cigarCount; ++i)
            {
                CigarOperator operator = cigarElements.get(i).getOperator();

                if(operator == lastOperator)
                    return false;

                if(lastOperator == S && operator != M)
                    return false;

                if(i == cigarCount - 1 && operator != M && operator != S)
                    return false;

                lastOperator = operator;
            }
        }
        else if(cigarElements.get(0).getOperator() != M)
        {
            return false;
        }

        return true;
    }
}
