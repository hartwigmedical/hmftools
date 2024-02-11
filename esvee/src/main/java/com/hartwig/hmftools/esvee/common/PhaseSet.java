package com.hartwig.hmftools.esvee.common;

import static java.lang.String.format;

import java.util.List;

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

    public void addAssemblyLink(final AssemblyLink link, int index)
    {
        mAssemblyLinks.add(index, link);
    }

    public List<AssemblyLink> assemblyLinks() { return mAssemblyLinks; }

    public String toString() { return format("id(%d) links(%d)", mId, mAssemblyLinks.size()); }
}
