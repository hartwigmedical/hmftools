package com.hartwig.hmftools.esvee.alignment;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;

public class AssemblyAlignment
{
    private final int mId;
    private final List<JunctionAssembly> mAssemblies;
    private final AssemblyLink mAssemblyLink;

    private final String mFullSequence;
    private final int mFullSequenceLength;

    private final List<Breakend> mBreakends;

    public AssemblyAlignment(final int id, final JunctionAssembly assembly) { this(id, assembly, null); }

    public AssemblyAlignment(final int id, final AssemblyLink assemblyLink) { this(id, null, assemblyLink); }

    private AssemblyAlignment(final int id, final JunctionAssembly assembly, final AssemblyLink assemblyLink)
    {
        mId = id;

        if(assembly != null)
        {
            mAssemblies = List.of(assembly);
            mAssemblyLink = null;
        }
        else
        {
            mAssemblies = List.of(assemblyLink.first(), assemblyLink.second());
            mAssemblyLink = assemblyLink;
        }

        mBreakends = Lists.newArrayList();

        mAssemblies.forEach(x -> x.setAssemblyAlignmentInfo(info()));

        mFullSequence = formFullSequence();
        mFullSequenceLength = mFullSequence.length();
    }

    public int id() { return mId; }

    public List<JunctionAssembly> assemblies() { return mAssemblies; }
    public JunctionAssembly first() { return mAssemblies.get(0); }
    public JunctionAssembly second() { return mAssemblies.size() > 1 ? mAssemblies.get(1) : null; }

    public List<Breakend> breakends() { return mBreakends; }
    public void addBreakend(final Breakend breakend) { mBreakends.add(breakend); }

    public String fullSequence() { return mFullSequence; }
    public int fullSequenceLength() { return mFullSequenceLength; }

    public void setOutcome(final AlignmentOutcome outcome)
    {
        mAssemblies.forEach(x -> x.setAlignmentOutcome(outcome));
    }

    public String assemblyIds()
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

    public int refBaseLength()
    {
        int refBaseLength = mAssemblies.get(0).refBaseLength();

        if(mAssemblies.size() > 1)
            refBaseLength += mAssemblies.get(1).refBaseLength();

        return refBaseLength;
    }

    private String formFullSequence()
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

        fullSequence.append(firstReversed ? Nucleotides.reverseComplementBases(first.formRefBaseSequence()) : first.formRefBaseSequence());

        if(!insertedBases.isEmpty())
            fullSequence.append(insertedBases);

        String secondSequence = secondReversed ? Nucleotides.reverseComplementBases(second.formRefBaseSequence()) : second.formRefBaseSequence();

        fullSequence.append(overlapLength > 0 ? secondSequence.substring(overlapLength) : secondSequence);

        return fullSequence.toString();
    }

    public String assemblyCigar()
    {
        if(mAssemblies.size() == 1)
        {
            JunctionAssembly assembly = mAssemblies.get(0);

            if(assembly.isForwardJunction())
                return format("%dM%dS", assembly.refBaseLength(), assembly.extensionLength());
            else
                return format("%dS%dM", assembly.extensionLength(), assembly.refBaseLength());
        }

        int insertedBases = mAssemblyLink.insertedBases().length();
        int overlapLength = mAssemblyLink.overlapBases().length();

        StringBuilder sb = new StringBuilder();
        sb.append(format("%dM", mAssemblyLink.first().refBaseLength()));

        if(insertedBases > 0)
            sb.append(format("%dI", insertedBases));

        int secondLength = mAssemblyLink.second().refBaseLength() - overlapLength;
        sb.append(format("%dM", secondLength));

        return sb.toString();
    }

    public String toString() { return info(); }
}
