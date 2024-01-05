package com.hartwig.hmftools.esvee.old;

import static com.hartwig.hmftools.esvee.read.ReadUtils.flipRead;

import com.hartwig.hmftools.esvee.read.Read;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.util.SequenceUtil;

public class ExtendedAssembly extends SupportedAssembly implements TrimmableAssembly<ExtendedAssembly>
{
    public final SupportedAssembly Source;

    public ExtendedAssembly(final String name, final String assembly, final SupportedAssembly source)
    {
        super(name, assembly);
        Source = source;
    }

    @Nullable
    @Override
    public ExtendedAssembly trim(final int removeLeft, final int removeRight)
    {
        final int newLength = getLength() - removeLeft - removeRight;
        if(newLength <= 0)
            return null;

        final String newBases = Assembly.substring(removeLeft, removeLeft + newLength);

        final ExtendedAssembly newAssembly = new ExtendedAssembly(Name, newBases, Source);

        for(ReadSupport support : readSupport())
        {
            final int newOffset = support.Index - removeLeft;
            if(newOffset >= newLength)
                continue;

            newAssembly.addEvidenceAt(support.Read, newOffset);
        }

        return newAssembly;
    }

    public ExtendedAssembly flipStrand()
    {
        final String assembly = SequenceUtil.reverseComplement(Assembly);
        final ExtendedAssembly flipped = new ExtendedAssembly(Name, assembly, Source);

        for(ReadSupport support : readSupport())
        {
            int initialReadLength = support.Read.getLength();
            Read flippedRead = flipRead(support.Read);

            flipped.addEvidenceAt(flippedRead,getLength() - support.Index - initialReadLength);
        }
        flipped.recalculateBaseQuality();
        return flipped;
    }
}
