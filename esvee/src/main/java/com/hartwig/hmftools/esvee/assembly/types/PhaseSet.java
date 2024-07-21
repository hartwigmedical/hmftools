package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.String.format;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.Orientation;

public class PhaseSet
{
    private int mId;
    private final List<AssemblyLink> mAssemblyLinks;
    private final List<JunctionAssembly> mAssemblies;

    public PhaseSet(final AssemblyLink link)
    {
        mId = -1;
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
            return;
        }

        int assemblyIndex = index == 0 ? 0 : mAssemblies.size();

        if(!mAssemblies.contains(link.first()))
        {
            mAssemblies.add(assemblyIndex, link.first());
        }

        if(!mAssemblies.contains(link.second()))
        {
            mAssemblies.add(assemblyIndex, link.second());
        }
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

    public Orientation assemblyOrientation(final JunctionAssembly assembly)
    {
        // the first assembly is defined as forward, meaning facing up the chain and each successive junction is alternating
        Orientation assemblyOrientation = Orientation.FORWARD;

        for(int i = 0; i < mAssemblies.size(); ++i)
        {
            if(mAssemblies.get(i) == assembly)
                return assemblyOrientation;

            assemblyOrientation = assemblyOrientation.opposite();
        }

        return null;
    }

    public boolean assembliesFaceInPhaseSet(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        Orientation assemblyOrientation = Orientation.FORWARD;

        int assemblyIndex1 = -1;
        Orientation assemblyOrientation1 = null;
        int assemblyIndex2 = -1;
        Orientation assemblyOrientation2 = null;

        for(int i = 0; i < mAssemblies.size(); ++i)
        {
            if(mAssemblies.get(i) == assembly1)
            {
                assemblyIndex1 = i;
                assemblyOrientation1 = assemblyOrientation;

                if(assemblyOrientation2 != null)
                    break;
            }
            else if(mAssemblies.get(i) == assembly2)
            {
                assemblyIndex2 = i;
                assemblyOrientation2 = assemblyOrientation;

                if(assemblyOrientation1 != null)
                    break;
            }

            assemblyOrientation = assemblyOrientation.opposite();
        }

        if(assemblyOrientation1 == null || assemblyOrientation2 == null)
            return false;

        if(assemblyOrientation1 == assemblyOrientation2)
            return false;

        // the read of the lower assembly (by index) faces up and vice versa
        if(assemblyIndex1 < assemblyIndex2)
        {
            return assemblyOrientation1.isForward();
        }
        else
        {
            return assemblyOrientation2.isForward();
        }
    }

    public String toString() { return format("id(%d) links(%d)", mId, mAssemblyLinks.size()); }
}
