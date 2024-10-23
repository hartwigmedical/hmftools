package com.hartwig.hmftools.esvee.assembly.phase;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.assembly.phase.ExtensionType.LOCAL_DEL_DUP;

import java.util.Comparator;

import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;

public class ExtensionCandidate
{
    public final ExtensionType Type;
    public final JunctionAssembly Assembly;
    public final JunctionAssembly SecondAssembly;
    public final AssemblyLink Link;

    public int SupportCount;

    public String ExtraInfo;
    public final Object Extender;

    private final boolean mIsValid;
    private boolean mSelected;

    public ExtensionCandidate(final ExtensionType type, final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        this(type, null, assembly1, assembly2);
    }

    public ExtensionCandidate(final ExtensionType type, final AssemblyLink assemblyLink)
    {
        this(type, assemblyLink, assemblyLink.first(), assemblyLink.second());
    }

    public ExtensionCandidate(
            final ExtensionType type, final AssemblyLink assemblyLink, final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        Type = type;
        Link = assemblyLink;
        mIsValid = assemblyLink != null;
        Assembly = assembly1;
        SecondAssembly = assembly2;
        SupportCount = 0;
        ExtraInfo = "";
        Extender = null;
        mSelected = false;
    }

    public ExtensionCandidate(final ExtensionType type, final JunctionAssembly assembly, final Object extender, final int supportCount)
    {
        Type = type;
        mIsValid = true;
        Assembly = assembly;
        Extender = extender;
        SupportCount = supportCount;
        ExtraInfo = "";

        Link = null;
        SecondAssembly = null;
        mSelected = false;
    }

    public boolean selected() { return mSelected; }
    public void markSelected() { mSelected = true; }

    public boolean matchesAssemblies(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        return (Assembly == assembly1 && SecondAssembly == assembly2) || (Assembly == assembly2 && SecondAssembly == assembly1);
    }

    public boolean isValid()
    {
        if(!mIsValid)
            return false;

        return SupportCount > 0 || Type == LOCAL_DEL_DUP;
    }

    protected static class StandardComparator implements Comparator<ExtensionCandidate>
    {
        @Override
        public int compare(final ExtensionCandidate first, final ExtensionCandidate second)
        {
            if(first.SupportCount != second.SupportCount)
                return first.SupportCount > second.SupportCount ? -1 : 1;

            if(first.Type == ExtensionType.SPLIT_LINK && second.Type == ExtensionType.SPLIT_LINK)
            {
                int firstJuncOffset = max(first.Link.insertedBases().length(), first.Link.overlapBases().length());
                int secondJuncOffset = max(first.Link.insertedBases().length(), first.Link.overlapBases().length());

                if(firstJuncOffset != secondJuncOffset)
                    return firstJuncOffset < secondJuncOffset ? -1 : 1; // favour a more precise link
            }

            // revert to type of link
            return Integer.compare(-first.Type.ordinal(), -second.Type.ordinal());
        }
    }

    protected static class LocalLinkComparator implements Comparator<ExtensionCandidate>
    {
        @Override
        public int compare(final ExtensionCandidate first, final ExtensionCandidate second)
        {
            if(first.Type != second.Type)
                return Integer.compare(-first.Type.ordinal(), -second.Type.ordinal());

            return Integer.compare(-first.SupportCount, -second.SupportCount);
        }
    }

    public String toString()
    {
        if(Link != null)
        {
            return format("%s %s link(%s) support(%d) %s",
                    mSelected ? "selected" : "candidate", Type, Link, SupportCount, ExtraInfo);
        }
        else
        {
            return format("%s %s assembly(%s) support(%d) %s",
                    mSelected ? "selected" : "candidate", Type, Assembly, SupportCount, ExtraInfo);
        }
    }
}
