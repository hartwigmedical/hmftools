package com.hartwig.hmftools.esvee.read;

import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;

public final class CigarUtils
{
    public static Cigar trimLeft(final Cigar cigar, final int count)
    {
        final List<CigarElement> elements = new ArrayList<>(cigar.getCigarElements());
        int toClip = count;
        while(toClip > 0)
        {
            final CigarElement element = elements.get(0);
            if(element.getLength() > toClip)
            {
                elements.set(0, new CigarElement(element.getLength() - toClip, element.getOperator()));
                break;
            }
            elements.remove(0);
            toClip -= element.getLength();
        }

        return new Cigar(elements);
    }

    public static Cigar trimRight(final Cigar cigar, final int count)
    {
        final List<CigarElement> elements = new ArrayList<>(cigar.getCigarElements());
        int toClip = count;
        while(toClip > 0)
        {
            final CigarElement element = elements.get(elements.size() - 1);
            if(element.getLength() > toClip)
            {
                elements.set(elements.size() - 1, new CigarElement(element.getLength() - toClip, element.getOperator()));
                break;
            }
            elements.remove(elements.size() - 1);
            toClip -= element.getLength();
        }

        return new Cigar(elements);
    }
}
