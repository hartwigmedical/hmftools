package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.flipOrientation;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

public class PhaseSet
{
    private int mId;
    private final List<AssemblyLink> mAssemblyLinks;
    private final List<JunctionAssembly> mAssemblies;
    private final List<Byte> mAssemblyOrientations;

    public PhaseSet(final AssemblyLink link)
    {
        mId = -1;
        mAssemblyOrientations = Lists.newArrayList();
        mAssemblyLinks = Lists.newArrayList();
        mAssemblies = Lists.newArrayList();
        addAssemblyLink(link, 0);
    }

    public void setId(int id) { mId = id; }
    public int id() { return mId; }

    public void addAssemblyLinkStart(final AssemblyLink link) { addAssemblyLink(link, 0); }
    public void addAssemblyLinkEnd(final AssemblyLink link) { addAssemblyLink(link, mAssemblyLinks.size()); }

    private void addAssemblyLink(final AssemblyLink link, int index)
    {
        mAssemblyLinks.add(index, link);

        if(mAssemblies.isEmpty())
        {
            mAssemblies.add(link.first());
            mAssemblies.add(link.second());
            mAssemblyOrientations.add(link.first().junction().Orientation);
            mAssemblyOrientations.add(link.second().junction().Orientation);
            return;
        }

        int assemblyIndex = index == 0 ? 0 : mAssemblies.size();

        if(!mAssemblies.contains(link.first()))
        {
            addAssembly(link.first(), assemblyIndex, link.type() == LinkType.FACING);
        }

        if(!mAssemblies.contains(link.second()))
        {
            addAssembly(link.second(), assemblyIndex, link.type() == LinkType.FACING);
        }
    }

    private void addAssembly(final JunctionAssembly assembly, int index, boolean isFacingLink)
    {
        byte linkOrientation;

        if(mAssemblies.isEmpty() || !isFacingLink)
        {
            linkOrientation = assembly.junction().Orientation;
        }
        else
        {
            // opposite to the facing assembly it links with
            if(index == 0)
                linkOrientation = flipOrientation(mAssemblies.get(0).junction().Orientation);
            else
                linkOrientation = flipOrientation(mAssemblies.get(mAssemblies.size() - 1).junction().Orientation);
        }

        mAssemblies.add(index, assembly);
        mAssemblyOrientations.add(index, linkOrientation);
    }

    public List<AssemblyLink> findAssemblyLinks(final JunctionAssembly assembly)
    {
        return mAssemblyLinks.stream().filter(x -> x.hasAssembly(assembly)).collect(Collectors.toList());
    }

    public AssemblyLink findSplitLink(final JunctionAssembly assembly)
    {
        return mAssemblyLinks.stream().filter(x -> x.hasAssembly(assembly)).filter(x -> x.type() == LinkType.SPLIT).findFirst().orElse(null);
    }

    public boolean hasAssembly(final JunctionAssembly assembly)
    {
        return mAssemblies.contains(assembly);
    }

    public List<AssemblyLink> assemblyLinks() { return mAssemblyLinks; }
    public List<JunctionAssembly> assemblies() { return mAssemblies; }
    public int linkCount() { return mAssemblyLinks.size(); }

    public int assemblyIndex(final JunctionAssembly assembly)
    {
        for(int i = 0; i < mAssemblies.size(); ++i)
        {
            if(mAssemblies.get(i) == assembly)
                return i;

        }

        return -1;
    }

    public byte assemblyOrientation(final JunctionAssembly assembly)
    {
        for(int i = 0; i < mAssemblies.size(); ++i)
        {
            if(mAssemblies.get(i) == assembly)
                return mAssemblyOrientations.get(i);

        }

        return -1;
    }

    public static boolean readsFaceInPhaseSet(
            final JunctionAssembly assembly1, int assemblyIndex1, byte assemblyOrientation1, final SupportRead read1,
            final JunctionAssembly assembly2, int assemblyIndex2, byte assemblyOrientation2, final SupportRead read2)
    {
        if(read1 == null || read2 == null)
            return false;

        byte adjustedReadOrientation1 = assembly1.junction().Orientation == assemblyOrientation1 ?
                read1.orientation() : flipOrientation(read1.orientation());

        byte adjustedReadOrientation2 = assembly2.junction().Orientation == assemblyOrientation2 ?
                read2.orientation() : flipOrientation(read2.orientation());

        if(adjustedReadOrientation1 == adjustedReadOrientation2)
            return false;

        if(assemblyIndex1 < assemblyIndex2)
        {
            return adjustedReadOrientation1 == POS_ORIENT;
        }
        else
        {
            return adjustedReadOrientation2 == POS_ORIENT;
        }
    }

    public String toString() { return format("id(%d) links(%d)", mId, mAssemblyLinks.size()); }
}
