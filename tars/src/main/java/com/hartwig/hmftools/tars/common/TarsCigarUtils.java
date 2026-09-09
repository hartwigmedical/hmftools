package com.hartwig.hmftools.tars.common;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public final class TarsCigarUtils
{
    private TarsCigarUtils()
    {
    }

    public static Cigar normalize(final Cigar cigar)
    {
        return new Cigar(normalize(cigar.getCigarElements()));
    }

    public static List<CigarElement> normalize(final List<CigarElement> elements)
    {
        return mergeAdjacentSameOp(dropZeroLength(elements));
    }

    public static List<CigarElement> dropZeroLength(final List<CigarElement> elements)
    {
        List<CigarElement> result = new ArrayList<>(elements.size());
        for(CigarElement element : elements)
        {
            if(element.getLength() > 0)
            {
                result.add(element);
            }
        }
        return result;
    }

    public static List<CigarElement> mergeAdjacentSameOp(final List<CigarElement> elements)
    {
        List<CigarElement> merged = new ArrayList<>(elements.size());
        for(CigarElement element : elements)
        {
            if(!merged.isEmpty() && merged.get(merged.size() - 1).getOperator() == element.getOperator())
            {
                CigarElement previous = merged.remove(merged.size() - 1);
                merged.add(new CigarElement(previous.getLength() + element.getLength(), element.getOperator()));
            }
            else
            {
                merged.add(element);
            }
        }
        return merged;
    }

    public static boolean isMatchOrEqualOp(final CigarOperator op)
    {
        return op == CigarOperator.M || op == CigarOperator.EQ;
    }

    public static int terminalMatchedRun(final List<CigarElement> elements, final boolean fromEnd)
    {
        int start = fromEnd ? elements.size() - 1 : 0;
        int step = fromEnd ? -1 : 1;
        for(int i = start; i >= 0 && i < elements.size(); i += step)
        {
            CigarOperator op = elements.get(i).getOperator();
            if(op == CigarOperator.S)
            {
                continue;
            }
            if(op.isAlignment())
            {
                return elements.get(i).getLength();
            }
            return 0;
        }
        return 0;
    }

    public static boolean indelAdjacentToTerminalSoftClip(final List<CigarElement> elements, final boolean leadingSide)
    {
        if(elements.size() < 2)
        {
            return false;
        }
        if(leadingSide)
        {
            if(elements.get(0).getOperator() != CigarOperator.S)
            {
                return false;
            }
            CigarOperator next = elements.get(1).getOperator();
            return next == CigarOperator.I || next == CigarOperator.D;
        }

        int last = elements.size() - 1;
        if(elements.get(last).getOperator() != CigarOperator.S)
        {
            return false;
        }
        CigarOperator previous = elements.get(last - 1).getOperator();
        return previous == CigarOperator.I || previous == CigarOperator.D;
    }

    public static Cigar clampLeadingReferenceToSoftClip(final Cigar cigar, final int overhang)
    {
        List<CigarElement> elements = cigar.getCigarElements();
        int existingLeadingSoftClip = leftSoftClipLength(cigar);
        int index = existingLeadingSoftClip > 0 ? 1 : 0;

        int remainingOverhang = overhang;
        int clippedReadBases = 0;
        List<CigarElement> tail = null;

        while(index < elements.size())
        {
            CigarElement element = elements.get(index);
            CigarOperator op = element.getOperator();

            if(!op.consumesReferenceBases())
            {
                if(op.consumesReadBases())
                {
                    clippedReadBases += element.getLength();
                }
                ++index;
                continue;
            }

            if(element.getLength() <= remainingOverhang)
            {
                remainingOverhang -= element.getLength();
                if(op.consumesReadBases())
                {
                    clippedReadBases += element.getLength();
                }
                ++index;

                if(remainingOverhang == 0)
                {
                    tail = new ArrayList<>(elements.subList(index, elements.size()));
                    break;
                }
                continue;
            }

            if(op.consumesReadBases())
            {
                clippedReadBases += remainingOverhang;
            }

            tail = new ArrayList<>();
            tail.add(new CigarElement(element.getLength() - remainingOverhang, op));
            tail.addAll(elements.subList(index + 1, elements.size()));
            remainingOverhang = 0;
            break;
        }

        if(remainingOverhang > 0 || tail == null || tail.isEmpty() || tail.get(0).getOperator() != CigarOperator.M)
        {
            return null;
        }

        List<CigarElement> out = new ArrayList<>(tail.size() + 1);
        out.add(new CigarElement(existingLeadingSoftClip + clippedReadBases, CigarOperator.S));
        out.addAll(tail);
        return new Cigar(out);
    }

    public static Cigar clampTrailingReferenceToSoftClip(final Cigar cigar, final int overhang)
    {
        List<CigarElement> reversed = new ArrayList<>(cigar.getCigarElements());
        Collections.reverse(reversed);

        Cigar clamped = clampLeadingReferenceToSoftClip(new Cigar(reversed), overhang);
        if(clamped == null)
        {
            return null;
        }

        List<CigarElement> out = new ArrayList<>(clamped.getCigarElements());
        Collections.reverse(out);
        return new Cigar(out);
    }
}
