package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.AssemblyConstants.PROXIMATE_JUNCTION_DISTANCE;
import static com.hartwig.hmftools.esvee.assembly.read.ReadFilters.recordSoftClipsAtJunction;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;

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
                || newAssembly.junction().Orient != assembly.junction().Orient)
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
                if(selectFirstAssembly(assembly, newAssembly))
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


    private static boolean selectFirstAssembly(final JunctionAssembly first, final JunctionAssembly second)
    {
        if(first.supportCount() != second.supportCount())
            return first.supportCount() > second.supportCount();

        // take the one with the most precise support
        long preciseFirst = first.support().stream().filter(x -> recordSoftClipsAtJunction(x.cachedRead(), first.junction())).count();
        long preciseSecond = second.support().stream().filter(x -> recordSoftClipsAtJunction(x.cachedRead(), second.junction())).count();

        return preciseFirst >= preciseSecond;
    }
}
