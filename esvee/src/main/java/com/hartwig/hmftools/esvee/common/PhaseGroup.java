package com.hartwig.hmftools.esvee.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;

public class PhaseGroup
{
    private int mId;
    private final List<JunctionAssembly> mAssemblies;
    private final List<DiscordantGroup> mDiscordantGroups;
    private final List<PhaseSet> mPhaseSets;
    private final List<AssemblyLink> mSecondarySplitLinks;

    public PhaseGroup(final JunctionAssembly first, @Nullable final JunctionAssembly second)
    {
        mId = -1;
        mAssemblies = Lists.newArrayList(first);
        mPhaseSets = Lists.newArrayList();
        mSecondarySplitLinks = Lists.newArrayList();
        mDiscordantGroups = Lists.newArrayList();

        first.setPhaseGroup(this, second != null ? second.junction().coords() : "");

        if(second != null)
        {
            mAssemblies.add(second);
            second.setPhaseGroup(this, first.junction().coords());
        }
    }

    public void setId(int id) { mId = id; }
    public int id() { return mId; }

    public List<JunctionAssembly> assemblies() { return mAssemblies; }
    public List<DiscordantGroup> discordantGroups() { return mDiscordantGroups; }
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

    public boolean isSolo() { return mAssemblies.size() == 1; }
    public int assemblyCount() { return mAssemblies.size(); }

    public void transferAssemblies(final PhaseGroup other)
    {
        for(JunctionAssembly assembly : other.assemblies())
        {
            mAssemblies.add(assembly);
            assembly.setPhaseGroup(this, "transfer");
        }
    }

    public void addAssembly(final JunctionAssembly assembly, final JunctionAssembly linkingAssembly)
    {
        mAssemblies.add(assembly);

        if(assembly.phaseGroup() != null)
        {
            SV_LOGGER.error("assembly({}) adding to additional phase group", assembly);
            System.exit(1);
        }

        assembly.setPhaseGroup(this, linkingAssembly != null ? linkingAssembly.junction().coords() : "");
    }

    public void addDiscordantGroup(final DiscordantGroup discordantGroup)
    {
        mDiscordantGroups.add(discordantGroup);

        /*
        if(discordantGroup.phaseGroup() != null)
        {
            // not sure yet whether this matters
            SV_LOGGER.warn("discGroup({}) adding to additional phase group", discordantGroup);
        }

        discordantGroup.setPhaseGroup(this);
        */
    }

    public String toString() { return format("id(%d) assemblies(%d) discGroups(%d)", mId, mAssemblies.size(), mDiscordantGroups.size()); }
}
