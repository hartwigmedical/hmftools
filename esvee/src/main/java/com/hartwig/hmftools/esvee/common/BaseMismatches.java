package com.hartwig.hmftools.esvee.common;

import static java.lang.String.format;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.esvee.read.Read;

public class BaseMismatches
{
    public final BaseMismatch[] Mismatches; // could consider a map or list since surely 1 value is most common, but still an overhead

    public BaseMismatches(final BaseMismatch baseMismatch)
    {
        Mismatches = new BaseMismatch[] { null, null, null, null};
        Mismatches[baseMismatch.Base.ordinal()] = baseMismatch;
    }

    public int mismatchCount() { return (int)Arrays.stream(Mismatches).filter(x -> x != null).count(); }

    public void addMismatch(final byte base, final Read read, final byte qual)
    {
        BaseType baseType = BaseType.from(base);

        if(Mismatches[baseType.ordinal()] == null)
        {
            Mismatches[baseType.ordinal()] = new BaseMismatch(base, read, qual);
        }
        else
        {
            Mismatches[baseType.ordinal()].QualTotal += qual;
            Mismatches[baseType.ordinal()].Reads.add(read);
        }
    }

    public String toString() { return format("mismatches(%d)", mismatchCount()); }

    @VisibleForTesting
    public List<BaseMismatch> baseMismatches() { return Arrays.stream(Mismatches).filter(x -> x != null).collect(Collectors.toList()); }
}
