package com.hartwig.hmftools.esvee.common;

import static java.lang.String.format;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.read.Read;

public class BaseMismatches
{
    public final BaseMismatch[] Mismatches;

    public BaseMismatches(final BaseMismatch baseMismatch)
    {
        Mismatches = new BaseMismatch[] { null, null, null, null};
        Mismatches[baseMismatch.Base.ordinal()] = baseMismatch;
    }

    public void addMismatch(final byte base, final Read read, final int qual)
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

    /*
    public void addMismatch(final BaseMismatch baseMismatch)
    {
        if(Mismatches[baseMismatch.Base.ordinal()] == null)
        {
            Mismatches[baseMismatch.Base.ordinal()] = baseMismatch;
        }
        else
        {
            Mismatches[baseMismatch.Base.ordinal()].QualTotal += baseMismatch.QualTotal;
        }
    }

    private BaseMismatch getOrCreate(final BaseMismatch baseMismatch)
    {
        if(Mismatches[baseMismatch.Base.ordinal()] == null)
        {
            Mismatches[baseMismatch.Base.ordinal()] = baseMismatch;
        }
        else
        {

        }

        return Mismatches[baseMismatch.Base.ordinal()];
    }
    */
}
