package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.assembly.read.Read;

public class BaseMismatch
{
    public final BaseType Base;
    public final List<Read> Reads;
    public byte MaxQual;
    public int QualTotal;

    public BaseMismatch(final byte base, final Read read, final byte qual)
    {
        Base = BaseType.from(base);
        Reads = Lists.newArrayList(read);
        MaxQual = qual;
        QualTotal = qual;
    }

    public byte base() { return Base.Byte; }

    public void addRead(final Read read, final byte qual)
    {
        Reads.add(read);
        MaxQual = maxQual(qual, MaxQual);
        QualTotal += qual;
    }

    public static byte maxQual(final byte first, final byte second) { return first > second ? first : second; }

    public String toString() { return format("%s reads(%d) totalQual(%d)", Base, Reads.size(), QualTotal); }
}
