package com.hartwig.hmftools.esvee.assembly.phase;

import static java.lang.String.format;

import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;

public class ExtensionCandidate
{
    public final ExtensionType Type;
    public final JunctionAssembly Assembly;
    public final JunctionAssembly SecondAssembly;
    public final AssemblyLink Link;

    public int AssemblyMatchedSupport;
    public int SecondAssemblyMatchedSupport;
    public int AssemblyCandidateReads;
    public int SecondAssemblyCandidateReads;

    public String ExtraInfo;
    public final Object Extender;

    public final boolean mIsValid;

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
        AssemblyMatchedSupport = 0;
        SecondAssemblyMatchedSupport = 0;
        AssemblyCandidateReads = 0;
        SecondAssemblyCandidateReads = 0;
        ExtraInfo = "";
        Extender = null;
    }

    public ExtensionCandidate(final ExtensionType type, final JunctionAssembly assembly, final Object extender, final int candidates)
    {
        Type = type;
        mIsValid = true;
        Assembly = assembly;
        Extender = extender;
        AssemblyCandidateReads = candidates;
        AssemblyMatchedSupport = 0;
        ExtraInfo = "";

        Link = null;
        SecondAssembly = null;
        SecondAssemblyMatchedSupport = 0;
        SecondAssemblyCandidateReads = 0;
    }

    public int totalSupport()
    {
        return AssemblyMatchedSupport + SecondAssemblyMatchedSupport + AssemblyCandidateReads + SecondAssemblyCandidateReads;
    }

    public boolean isValid() { return mIsValid && totalSupport() > 0; }

    public String toString()
    {
        if(Link != null)
        {
            return format("%s link(%s) support first(s=%d c=%d) second(s=%d c=%d) %s",
                    Type, Link, AssemblyMatchedSupport, AssemblyCandidateReads, SecondAssemblyMatchedSupport,
                    SecondAssemblyCandidateReads, ExtraInfo);
        }
        else
        {
            return format("%s assembly(%s) support(s=%d c=%d) %s",
                    Type, Assembly, AssemblyMatchedSupport, AssemblyCandidateReads, ExtraInfo);
        }
    }
}
