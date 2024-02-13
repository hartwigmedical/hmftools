package com.hartwig.hmftools.esvee.common;

import static java.lang.String.format;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

public class PhaseSet
{
    private int mId;
    private final List<AssemblyLink> mAssemblyLinks;

    public PhaseSet(final AssemblyLink link)
    {
        mId = -1;
        mAssemblyLinks = Lists.newArrayList(link);
    }

    public void setId(int id) { mId = id; }
    public int id() { return mId; }

    public int linkCount() { return mAssemblyLinks.size(); }

    public void addAssemblyStart(final AssemblyLink link) { addAssemblyLink(link, 0); }
    public void addAssemblyEnd(final AssemblyLink link) { addAssemblyLink(link, mAssemblyLinks.size()); }

    private void addAssemblyLink(final AssemblyLink link, int index)
    {
        mAssemblyLinks.add(index, link);
    }

    public List<AssemblyLink> findAssemblyLinks(final JunctionAssembly assembly)
    {
        return mAssemblyLinks.stream().filter(x -> x.hasAssembly(assembly)).collect(Collectors.toList());
    }

    public boolean hasAssembly(final JunctionAssembly assembly)
    {
        return mAssemblyLinks.stream().anyMatch(x -> x.hasAssembly(assembly));
    }

    public boolean hasMatchingAssemblyLink(final AssemblyLink link)
    {
        return mAssemblyLinks.stream().anyMatch(x -> x.matches(link));
    }

    public boolean hasMatchingAssembly(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        return mAssemblyLinks.stream().anyMatch(x -> x.matches(assembly1, assembly2));
    }

    public List<AssemblyLink> assemblyLinks() { return mAssemblyLinks; }

    public String toString() { return format("id(%d) links(%d)", mId, mAssemblyLinks.size()); }
}
