package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MERGE_READ_OVERLAP;
import static com.hartwig.hmftools.esvee.SvConstants.PROXIMATE_JUNCTION_DISTANCE;
import static com.hartwig.hmftools.esvee.SvConstants.PROXIMATE_JUNCTION_OVERLAP;

import java.util.List;

import com.google.common.collect.Lists;
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

    public static void dedupProximateAssemblies(final List<JunctionAssembly> existingAssemblies, final List<JunctionAssembly> newAssemblies)
    {
        if(newAssemblies.isEmpty() || existingAssemblies.isEmpty())
            return;

        // start with the most recent previous assemblies since they are added in order
        int index = existingAssemblies.size() - 1;
        int minPosition = newAssemblies.get(0).junction().Position - PROXIMATE_JUNCTION_DISTANCE;

        List<JunctionAssembly> removeExisting = Lists.newArrayList();

        while(index >= 0)
        {
            JunctionAssembly assembly = existingAssemblies.get(index);

            if(assembly.junction().Position < minPosition)
                break;

            int newIndex = 0;

            while(newIndex < newAssemblies.size())
            {
                JunctionAssembly newAssembly = newAssemblies.get(newIndex);

                if(newAssembly.junction() == assembly.junction()
                || newAssembly.junction().Orientation != assembly.junction().Orientation)
                {
                    ++newIndex;
                    continue;
                }

                if(!haveOverlapDistance(assembly, newAssembly) || !haveOverlappingReads(assembly, newAssembly))
                {
                    ++newIndex;
                    continue;
                }

                if(!SequenceCompare.matchedAssemblySequences(assembly, newAssembly))
                {
                    ++newIndex;
                    continue;
                }

                // take the assembly with the most read support
                if(assembly.supportCount() >= newAssembly.supportCount())
                {
                    assembly.addMergedAssembly();
                    newAssemblies.remove(newIndex);
                }
                else
                {
                    // keep this assembly and mark the existing one for removal
                    newAssembly.addMergedAssembly();
                    removeExisting.add(assembly);
                    ++newIndex;
                }
            }

            --index;
        }

        removeExisting.forEach(x -> existingAssemblies.remove(x));
    }

    private static boolean haveOverlapDistance(final JunctionAssembly lower, final JunctionAssembly upper)
    {
        int overlapDistance = min(upper.maxAlignedPosition(), lower.maxAlignedPosition())
                - max(upper.minAlignedPosition(), lower.minAlignedPosition());

        return overlapDistance >= PROXIMATE_JUNCTION_OVERLAP;
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
