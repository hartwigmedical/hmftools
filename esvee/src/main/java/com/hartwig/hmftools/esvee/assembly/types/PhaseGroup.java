package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.util.List;
import java.util.stream.Collectors;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;

public class PhaseGroup
{
    private int mId;
    private final List<JunctionAssembly> mAssemblies;
    private final List<JunctionAssembly> mDerivedAssemblies;
    private final List<PhaseSet> mPhaseSets;
    private final List<AssemblyLink> mSecondarySplitLinks;

    public PhaseGroup(final JunctionAssembly first, @Nullable final JunctionAssembly second)
    {
        mId = -1;
        mAssemblies = Lists.newArrayList(first);
        mDerivedAssemblies = Lists.newArrayList();
        mPhaseSets = Lists.newArrayList();
        mSecondarySplitLinks = Lists.newArrayList();

        first.setPhaseGroup(this);

        if(second != null)
        {
            mAssemblies.add(second);
            second.setPhaseGroup(this);
        }
    }

    public void setId(int id) { mId = id; }
    public int id() { return mId; }

    public List<JunctionAssembly> assemblies() { return mAssemblies; }
    public List<JunctionAssembly> derivedAssemblies() { return mDerivedAssemblies; }
    public List<PhaseSet> phaseSets() { return mPhaseSets; }

    public List<AssemblyLink> secondaryLinks() { return mSecondarySplitLinks; }

    public PhaseSet findPhaseSet(final JunctionAssembly assembly)
    {
        return mPhaseSets.stream().filter(x -> x.hasAssembly(assembly)).findFirst().orElse(null);
    }

    public List<AssemblyLink> findSecondarySplitLinks(final JunctionAssembly assembly)
    {
        return mSecondarySplitLinks.stream().filter(x -> x.hasAssembly(assembly)).collect(Collectors.toList());
    }

    public int assemblyCount() { return mAssemblies.size(); }

    public void transferAssemblies(final PhaseGroup other)
    {
        for(JunctionAssembly assembly : other.assemblies())
        {
            if(mAssemblies.contains(assembly))
            {
                SV_LOGGER.error("assembly({}) transferred from pg({})but already in phase group({})",
                        assembly, other, this);
                System.exit(1);
            }

            mAssemblies.add(assembly);

            assembly.setPhaseGroup(this);
        }
    }

    public void addAssembly(final JunctionAssembly assembly)
    {
        mAssemblies.add(assembly);

        if(assembly.phaseGroup() != null)
        {
            SV_LOGGER.error("assembly({}) adding to additional phase group", assembly);
            System.exit(1);
        }

        assembly.setPhaseGroup(this);
    }

    public void addDerivedAssembly(final JunctionAssembly assembly)
    {
        mDerivedAssemblies.add(assembly);
        addAssembly(assembly);
    }

    public String toString() { return format("id(%d) assemblies(%d)", mId, mAssemblies.size()); }
}
