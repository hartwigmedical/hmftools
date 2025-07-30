package com.hartwig.hmftools.redux.jitter;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.util.StringUtil;

public class RefGenomeMicrosatellite
{
    public static final Logger sLogger = LogManager.getLogger(RefGenomeMicrosatellite.class);

    public final ChrBaseRegion genomeRegion;
    public final byte[] unit;
    public final int numRepeat;
    public double mappability = Double.NaN;

    public RefGenomeMicrosatellite(final ChrBaseRegion genomeRegion, final byte[] unit)
    {
        this.genomeRegion = genomeRegion;
        this.unit = unit;
        this.numRepeat = genomeRegion.baseLength() / unit.length;
    }

    public RefGenomeMicrosatellite(final ChrBaseRegion genomeRegion, byte unit)
    {
        this(genomeRegion, new byte[] { unit });
    }

    public RefGenomeMicrosatellite(final String chromosome, int start, int end, final byte[] unit)
    {
        this(new ChrBaseRegion(chromosome, start, end), unit);
    }

    public RefGenomeMicrosatellite(final String chromosome, int start, int end, final byte unit)
    {
        this(new ChrBaseRegion(chromosome, start, end), unit);
    }

    public String chromosome()
    {
        return genomeRegion.chromosome();
    }

    public int referenceStart()
    {
        return genomeRegion.start();
    }

    public int referenceEnd()
    {
        return genomeRegion.end();
    }

    public int baseLength()
    {
        return genomeRegion.baseLength();
    }

    public String unitString()
    {
        return StringUtil.bytesToString(unit);
    }

    @Override
    public String toString()
    {
        return genomeRegion.toString() + ' ' + numRepeat + " x " + unitString();
    }
}
