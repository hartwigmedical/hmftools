package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.String.format;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

public class PhaseSet
{
    private int mId;
    private final List<AssemblyLink> mAssemblyLinks;
    private final List<JunctionAssembly> mAssemblies;

    public PhaseSet(final AssemblyLink link)
    {
        mId = -1;
        mAssemblyLinks = Lists.newArrayList(link);
        mAssemblies = Lists.newArrayList(link.first(), link.second());
    }

    public void setId(int id) { mId = id; }
    public int id() { return mId; }

    public int linkCount() { return mAssemblyLinks.size(); }

    public void addAssemblyStart(final AssemblyLink link) { addAssemblyLink(link, 0); }
    public void addAssemblyEnd(final AssemblyLink link) { addAssemblyLink(link, mAssemblyLinks.size()); }

    private void addAssemblyLink(final AssemblyLink link, int index)
    {
        mAssemblyLinks.add(index, link);

        if(!mAssemblies.contains(link.first()))
            mAssemblies.add(link.first());

        if(!mAssemblies.contains(link.second()))
            mAssemblies.add(link.second());
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

    public String toString() { return format("id(%d) links(%d)", mId, mAssemblyLinks.size()); }
}
