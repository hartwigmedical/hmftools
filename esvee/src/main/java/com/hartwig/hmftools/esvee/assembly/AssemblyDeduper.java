package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MERGE_READ_OVERLAP;

import java.util.List;

import com.hartwig.hmftools.esvee.common.AssemblySupport;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;

public class AssemblyDeduper
{
    public static void dedupJunctionAssemblies(final List<JunctionAssembly> assemblies)
    {
        if(assemblies.size() < 2)
            return;

        int index = 0;

        while(index < assemblies.size() - 1)
        {
            JunctionAssembly first = assemblies.get(index);

            int nextIndex = index + 1;

            while(nextIndex < assemblies.size())
            {
                JunctionAssembly second = assemblies.get(nextIndex);

                if(!haveOverlappingReads(first, second))
                {
                    ++nextIndex;
                    continue;
                }

                if(!SequenceCompare.matchedAssemblySequences(first, second))
                {
                    ++nextIndex;
                    continue;
                }

                // take the assembly with the most read support
                if(first.supportCount() >= second.supportCount())
                {
                    first.addMergedAssembly();
                }
                else
                {
                    second.addMergedAssembly();
                    assemblies.set(index, second);
                    first = second;
                }

                assemblies.remove(nextIndex);
            }

            ++index;
        }
    }

    private static boolean haveOverlappingReads(final JunctionAssembly first, final JunctionAssembly second)
    {
        int matchedReads = 0;

        for(AssemblySupport support : first.support())
        {
            if(second.support().stream().anyMatch(x -> x.read() == support.read()))
            {
                ++matchedReads;

                if(matchedReads >= PRIMARY_ASSEMBLY_MERGE_READ_OVERLAP)
                    return true;
            }
        }

        return false;
    }
}
