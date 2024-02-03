package com.hartwig.hmftools.esvee.common;

import java.util.List;

import com.google.common.collect.Lists;

public class PrimaryPhaseGroup
{
    private int mId;
    private final List<JunctionAssembly> mAssemblies;

    public PrimaryPhaseGroup(final JunctionAssembly first, final JunctionAssembly second)
    {
        mId = -1;
        mAssemblies = Lists.newArrayList(first, second);
        first.setPrimaryPhaseGroup(this);
        second.setPrimaryPhaseGroup(this);
    }

    public void setId(int id) { mId = id; }
    public int id() { return mId; }

    public int assemblyCount() { return mAssemblies.size(); }

    public void addAssembly(final JunctionAssembly assembly)
    {
        mAssemblies.add(assembly);
        assembly.setPrimaryPhaseGroup(this);
    }

    public List<JunctionAssembly> assemblies() { return mAssemblies; }
}
