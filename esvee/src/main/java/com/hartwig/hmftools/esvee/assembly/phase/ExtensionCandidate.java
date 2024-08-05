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

    public final String ExtraInfo;

    public ExtensionCandidate(
            final ExtensionType type, final AssemblyLink assemblyLink,  final int firstMatchedSupport, final int firstCandidateReads,
            final int secondMatchedSupport, final int secondCandidateReads, final String extraInfo)
    {
        Type = type;
        Link = assemblyLink;
        Assembly = assemblyLink.first();
        SecondAssembly = assemblyLink.second();
        AssemblyMatchedSupport = firstMatchedSupport;
        SecondAssemblyMatchedSupport = secondMatchedSupport;
        AssemblyCandidateReads = firstCandidateReads;
        SecondAssemblyCandidateReads = secondCandidateReads;
        ExtraInfo = extraInfo;
    }

    public ExtensionCandidate(final ExtensionType type, final JunctionAssembly assembly, final String extraInfo, final int candidates)
    {
        Type = type;
        Assembly = assembly;
        AssemblyCandidateReads = candidates;
        AssemblyMatchedSupport = 0;
        ExtraInfo = extraInfo;

        Link = null;
        SecondAssembly = null;
        SecondAssemblyMatchedSupport = 0;
        SecondAssemblyCandidateReads = 0;
    }

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
