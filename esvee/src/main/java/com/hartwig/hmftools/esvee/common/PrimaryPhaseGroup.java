package com.hartwig.hmftools.esvee.common;

import java.util.List;

import com.google.common.collect.Lists;

public class PrimaryPhaseGroup
{
    private final List<JunctionAssembly> mAssemblies;

    public PrimaryPhaseGroup(final JunctionAssembly first, final JunctionAssembly second)
    {
        mAssemblies = Lists.newArrayList(first, second);
        first.setPrimaryPhaseGroup(this);
        second.setPrimaryPhaseGroup(this);
    }

    public int assemblyCount() { return mAssemblies.size(); }

    public void addAssembly(final JunctionAssembly assembly)
    {
        mAssemblies.add(assembly);
        assembly.setPrimaryPhaseGroup(this);
    }

    public List<JunctionAssembly> assemblies() { return mAssemblies; }
}
