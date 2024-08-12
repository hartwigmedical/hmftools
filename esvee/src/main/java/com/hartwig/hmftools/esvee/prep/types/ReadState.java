package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.Math.max;
import static java.lang.String.format;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class ReadState
{
    public final int AlignedBaseLength;
    public final int SoftClipLengthLeft;
    public final int SoftClipLengthRight;
    public final int MaxIndelLength;

    private int mFilters;

    public ReadState(final SAMRecord record)
    {
        mFilters = 0;

        int alignedBaseLength = 0;
        int softClipLengthLeft = 0;
        int softClipLengthRight = 0;
        int maxIndelLength = 0;

        for(int i = 0; i < record.getCigar().getCigarElements().size(); ++i)
        {
            CigarElement element = record.getCigar().getCigarElements().get(i);

            switch(element.getOperator())
            {
                case M:
                    alignedBaseLength += element.getLength();
                    break;

                case S:
                    if(i == 0)
                        softClipLengthLeft = element.getLength();
                    else
                        softClipLengthRight = element.getLength();
                    break;

                case D:
                case I:
                    maxIndelLength = max(element.getLength(), maxIndelLength);
            }
        }

        AlignedBaseLength = alignedBaseLength;
        SoftClipLengthLeft = softClipLengthLeft;
        SoftClipLengthRight = softClipLengthRight;
        MaxIndelLength = maxIndelLength;
    }

    public void addFilter(final ReadFilterType filterType) { mFilters |= filterType.flag(); }
    public boolean hasFilter(final ReadFilterType filterType) { return (mFilters & filterType.flag()) != 0; }
    public void removefilter(final ReadFilterType filterType) { mFilters &= ~filterType.flag(); }
    public boolean unfiltered() { return mFilters == 0; }

    public int filters() { return mFilters; }

    public String toString()
    {
        return format("aligned(%d) softClip(%d - %d) maxIndel(%d) filters(%d)",
                AlignedBaseLength, SoftClipLengthLeft, SoftClipLengthRight, MaxIndelLength, mFilters);
    }
}
