package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.common.CommonUtils.formSvType;

import com.hartwig.hmftools.common.sv.StructuralVariantType;

public class AssemblyLink
{
    private final JunctionAssembly mFirst;
    private final JunctionAssembly mSecond;
    private final LinkType mType;

    private final String mInsertedBases;
    private final String mOverlapBases;

    public AssemblyLink(
            final JunctionAssembly first, final JunctionAssembly second, final LinkType type,
            final String insertedBases, final String overlapBases)
    {
        mFirst = first;
        mSecond = second;
        mType = type;
        mInsertedBases = insertedBases;
        mOverlapBases = overlapBases;
    }

    public LinkType type() { return mType; }
    public JunctionAssembly first() { return mFirst; }
    public JunctionAssembly second() { return mSecond; }

    public JunctionAssembly otherAssembly(final JunctionAssembly assembly) { return mFirst == assembly ? mSecond : mFirst; }

    public String insertedBases() { return mInsertedBases; }
    public String overlapBases() { return mOverlapBases; }

    public StructuralVariantType svType()
    {
        return formSvType(
                mFirst.junction().Chromosome, mSecond.junction().Chromosome, mFirst.junction().Position, mSecond.junction().Position,
                mFirst.junction().Orient, mSecond.junction().Orient, !mInsertedBases.isEmpty());
    }

    public boolean matches(final AssemblyLink other) { return matches(other.first(), other.second()); }

    public boolean matches(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        return (mFirst.equals(assembly1) && mSecond.equals(assembly2)) || (mFirst.equals(assembly2) && mSecond.equals(assembly1));
    }

    public boolean hasAssembly(final JunctionAssembly assembly) { return mFirst.equals(assembly) || mSecond.equals(assembly); }

    public int length()
    {
        if(mType == LinkType.FACING)
            return abs(mFirst.junction().Position - mSecond.junction().Position) + 1;

        return mFirst.junction().Chromosome.equals(mSecond.junction().Chromosome) ?
                abs(mFirst.junction().Position - mSecond.junction().Position) : 0;
    }

    public String toString()
    {
        return format("%s: %s - %s len(%d) extras(overlap=%d insert=%d)",
                mType, mFirst.junction().coords(), mSecond.junction().coords(), length(),
                mOverlapBases.length(), mInsertedBases.length());
    }
}
