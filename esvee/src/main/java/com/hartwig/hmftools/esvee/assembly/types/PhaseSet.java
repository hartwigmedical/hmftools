package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.SECONDARY;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.alignment.AssemblyAlignment;

public class PhaseSet
{
    private int mId;
    private final List<AssemblyLink> mAssemblyLinks;
    private final List<AssemblyLink> mSecondaryLinks;
    private final List<JunctionAssembly> mAssemblies;

    private AssemblyAlignment mAssemblyAlignment;
    private final List<PhaseSet> mMergedPhaseSets;
    private Integer mMergedPhaseSetId;

    public PhaseSet(final AssemblyLink link)
    {
        mId = -1;
        mAssemblyLinks = Lists.newArrayList();
        mAssemblies = Lists.newArrayList();
        mSecondaryLinks = Lists.newArrayList();
        addAssemblyLink(link, 0);

        mAssemblyAlignment = null;
        mMergedPhaseSets = Lists.newArrayList();
        mMergedPhaseSetId = null;
    }

    public void setId(int id) { mId = id; }
    public int id() { return mId; }

    public void addAssemblyLinkStart(final AssemblyLink link) { addAssemblyLink(link, 0); }
    public void addAssemblyLinkEnd(final AssemblyLink link) { addAssemblyLink(link, mAssemblyLinks.size()); }
    public void addSecondaryLink(final AssemblyLink link) { mSecondaryLinks.add(link); }

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

    public List<JunctionAssembly> assemblies() { return mAssemblies; }

    public List<AssemblyLink> assemblyLinks() { return mAssemblyLinks; }
    public List<AssemblyLink> secondaryLinks() { return mSecondaryLinks; }

    public boolean isSecondaryLineLink() { return mAssemblyLinks.size() == 1 && mAssemblies.stream().anyMatch(x -> x.outcome() == SECONDARY); }

    public AssemblyAlignment assemblyAlignment() { return mAssemblyAlignment; }
    public void setAssemblyAlignment(final AssemblyAlignment assemblyAlignment) { mAssemblyAlignment = assemblyAlignment; }

    public List<PhaseSet> mergedPhaseSets() { return mMergedPhaseSets; }
    public void mergePhaseSet(final PhaseSet phaseSet) { mMergedPhaseSets.add(phaseSet); }

    public boolean merged() { return mMergedPhaseSetId != null; }
    public int mergedPhaseSetId() { return mMergedPhaseSetId != null ? mMergedPhaseSetId : -1; }
    public void setMergedPhaseSetId(int phaseSetId) { mMergedPhaseSetId = phaseSetId; }

    public boolean hasAssembly(final JunctionAssembly assembly) { return mAssemblies.contains(assembly); }

    public List<AssemblyLink> findAssemblyLinks(final JunctionAssembly assembly)
    {
        return mAssemblyLinks.stream().filter(x -> x.hasAssembly(assembly)).collect(Collectors.toList());
    }

    public AssemblyLink findSplitLink(final JunctionAssembly assembly)
    {
        return mAssemblyLinks.stream().filter(x -> x.hasAssembly(assembly)).filter(x -> x.type() == LinkType.SPLIT).findFirst().orElse(null);
    }

    public boolean hasFacingLinks() { return mAssemblyLinks.stream().anyMatch(x -> x.type() == LinkType.FACING); }

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
