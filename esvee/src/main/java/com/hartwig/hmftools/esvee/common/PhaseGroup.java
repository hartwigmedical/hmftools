package com.hartwig.hmftools.esvee.common;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;

public class PhaseGroup
{
    private int mId;
    private final List<JunctionAssembly> mAssemblies;

    public PhaseGroup(final JunctionAssembly first, final JunctionAssembly second)
    {
        mId = -1;
        mAssemblies = Lists.newArrayList(first, second);
        first.setPhaseGroup(this);
        second.setPhaseGroup(this);
    }

    public PhaseGroup(final JunctionAssembly assembly)
    {
        mId = -1;
        mAssemblies = List.of(assembly);
    }

    public void setId(int id) { mId = id; }
    public int id() { return mId; }

    public boolean isSolo() { return mAssemblies.size() == 1; }
    public int assemblyCount() { return mAssemblies.size(); }

    public void addAssembly(final JunctionAssembly assembly)
    {
        mAssemblies.add(assembly);
        assembly.setPhaseGroup(this);
    }

    public List<JunctionAssembly> assemblies() { return mAssemblies; }

    public String toString() { return format("id(%d) assemblies(%d)", mId, mAssemblies.size()); }
}
