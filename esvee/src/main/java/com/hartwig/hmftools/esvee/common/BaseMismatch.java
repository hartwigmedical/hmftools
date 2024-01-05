package com.hartwig.hmftools.esvee.common;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.read.Read;

public class BaseMismatch
{
    public final BaseType Base;
    public final List<Read> Reads;
    public int QualTotal;

    public BaseMismatch(final byte base, final Read read, final int qual)
    {
        Base = BaseType.from(base);
        Reads = Lists.newArrayList(read);
        QualTotal = qual;
    }

    public byte base() { return Base.Byte; }

    public void addRead(final Read read, final int qual)
    {
        Reads.add(read);
        QualTotal += qual;
    }

    public String toString() { return format("%s reads(%d) totalQual(%d)", Base, Reads.size(), QualTotal); }
}
