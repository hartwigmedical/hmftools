package com.hartwig.hmftools.esvee.alignment;

import static java.lang.String.format;

import java.util.List;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.JunctionSequence;

public class AssemblyAlignment
{
    private final List<JunctionAssembly> mAssemblies;
    private final AssemblyLink mAssemblyLink;

    public AssemblyAlignment(final JunctionAssembly assembly)
    {
        mAssemblies = List.of(assembly);
        mAssemblyLink = null;
    }

    public AssemblyAlignment(final AssemblyLink assemblyLink)
    {
        mAssemblies = List.of(assemblyLink.first(), assemblyLink.second());
        mAssemblyLink = assemblyLink;
    }

    public String ids()
    {
        if(mAssemblies.size() == 1)
            return String.valueOf(mAssemblies.get(0).id());

        return format("%d;%d", mAssemblyLink.first().id(), mAssemblyLink.second().id());
    }

    public String info()
    {
        if(mAssemblies.size() == 1)
            return mAssemblies.get(0).junction().coords();

        return format("%s_%s", mAssemblyLink.first().junction().coords(), mAssemblyLink.second().junction().coords());
    }

    public StructuralVariantType svType()
    {
        return mAssemblies.size() == 1 ? StructuralVariantType.SGL : mAssemblyLink.svType();
    }

    public int svLength()
    {
        return mAssemblies.size() == 1 ? 0 : mAssemblyLink.length();
    }

    public String fullSequence()
    {
        if(mAssemblies.size() == 1)
            return mAssemblies.get(0).formFullSequence();

        // factor in orientation and overlapping or inserted bases
        JunctionAssembly first, second;
        boolean firstReversed = false;
        boolean secondReversed = false;

        if(mAssemblyLink.first().junction().Orientation != mAssemblyLink.second().junction().Orientation)
        {
            if(mAssemblyLink.first().junction().isForward())
            {
                first = mAssemblyLink.first();
                second = mAssemblyLink.second();

            }
            else
            {
                first = mAssemblyLink.second();
                second = mAssemblyLink.second();
            }
        }
        else
        {
            first = mAssemblyLink.first();
            second = mAssemblyLink.second();

            if(mAssemblyLink.first().junction().isForward())
                secondReversed = true;
            else
                firstReversed = true;
        }

        // now build a sequence of first followed by second
        int overlapLength = mAssemblyLink.overlapBases().length();
        String insertedBases = mAssemblyLink.insertedBases();

        int fullSequenceLength = first.baseLength() + second.baseLength() + insertedBases.length() - overlapLength;

        StringBuilder fullSequence = new StringBuilder(fullSequenceLength);

        fullSequence.append(firstReversed ? Nucleotides.reverseComplementBases(first.formFullSequence()) : first.formFullSequence());

        if(!insertedBases.isEmpty())
            fullSequence.append(insertedBases);

        String secondSequence = secondReversed ? Nucleotides.reverseComplementBases(second.formFullSequence()) : second.formFullSequence();

        fullSequence.append(overlapLength > 0 ? secondSequence.substring(overlapLength) : secondSequence);

        return fullSequence.toString();
    }

    public String toString() { return info(); }
}
