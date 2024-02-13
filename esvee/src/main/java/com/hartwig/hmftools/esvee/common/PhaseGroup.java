package com.hartwig.hmftools.esvee.common;

import static java.lang.String.format;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

public class PhaseGroup
{
    private int mId;
    private final List<JunctionAssembly> mAssemblies;
    private final List<PhaseSet> mPhaseSets;

    public PhaseGroup(final JunctionAssembly first, final JunctionAssembly second)
    {
        mId = -1;
        mAssemblies = Lists.newArrayList(first, second);
        mPhaseSets = Lists.newArrayList();
        first.setPhaseGroup(this);
        second.setPhaseGroup(this);
    }

    public PhaseGroup(final JunctionAssembly assembly)
    {
        mId = -1;
        mAssemblies = List.of(assembly);
        mPhaseSets = Collections.emptyList();
    }

    public void setId(int id)
    {
        mId = id;

        // also set phase set IDs
        int phaseSetId = 0;
        for(PhaseSet phaseSet : mPhaseSets)
        {
            phaseSet.setId(phaseSetId++);
        }
    }

    public int id() { return mId; }

    public List<PhaseSet> phaseSets() { return mPhaseSets; }

    public PhaseSet findPhaseSet(final JunctionAssembly assembly)
    {
        return mPhaseSets.stream().filter(x -> x.hasAssembly(assembly)).findFirst().orElse(null);
    }

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
