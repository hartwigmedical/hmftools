package com.hartwig.hmftools.redux.jitter;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.util.StringUtil;

public class RefGenomeMicrosatellite
{
    public final ChrBaseRegion Region;
    public final byte[] Unit;
    public final int RepeatCount;
    public double mMappability;

    public RefGenomeMicrosatellite(final ChrBaseRegion region, final byte[] unit)
    {
        Region = region;
        Unit = unit;
        RepeatCount = region.baseLength() / unit.length;
        mMappability = Double.NaN;
    }

    public RefGenomeMicrosatellite(final String chromosome, int start, int end, final byte[] unit)
    {
        this(new ChrBaseRegion(chromosome, start, end), unit);
    }

    public String chromosome()
    {
        return Region.chromosome();
    }
    public int referenceStart()
    {
        return Region.start();
    }
    public int referenceEnd()
    {
        return Region.end();
    }

    public double mappability() { return mMappability; }
    public void setMappability(double mappability) { mMappability = mappability; }

    public String unitString()
    {
        return StringUtil.bytesToString(Unit);
    }

    @Override
    public String toString()
    {
        return Region.toString() + ' ' + RepeatCount + "x" + unitString();
    }
}
