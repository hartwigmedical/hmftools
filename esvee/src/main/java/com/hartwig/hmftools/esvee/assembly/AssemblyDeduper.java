package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_DEDUP_JITTER_MAX_DIST;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PROXIMATE_JUNCTION_DISTANCE;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.recordSoftClipsAtJunction;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.common.IndelCoords;

public class AssemblyDeduper
{
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

                if(!canDedupAssemblies(assembly, newAssembly))
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

    private static boolean canDedupAssemblies(final JunctionAssembly first, final JunctionAssembly second)
    {
        if(!SequenceCompare.matchedAssemblySequences(first, second))
            return false;

        // check that one assembly doesn't support the second one existing if sufficient far apart
        int junctionDistance = abs(first.junction().Position - second.junction().Position);

        if(junctionDistance <= ASSEMBLY_DEDUP_JITTER_MAX_DIST)
            return true;

        IndelCoords firstIndelCoords = first.indelCoords();
        IndelCoords secondIndelCoords = second.indelCoords();

        for(SupportRead firstRead : first.support())
        {
            if(!firstRead.type().isSplitSupport())
                continue;

            if(secondIndelCoords != null && firstRead.indelCoords() != null && firstRead.indelCoords().matches(secondIndelCoords))
                return false;

            for(SupportRead secondRead : second.support())
            {
                if(secondRead.cachedRead() == firstRead.cachedRead())
                    return false;

                if(firstIndelCoords != null && secondRead.indelCoords() != null &&secondRead.indelCoords().matches(firstIndelCoords))
                    return false;
            }
        }

        return true;
    }

    private static boolean selectFirstAssembly(final JunctionAssembly first, final JunctionAssembly second)
    {
        if(first.indel() != second.indel())
            return first.indel();

        if(first.supportCount() != second.supportCount())
            return first.supportCount() > second.supportCount();

        // take the one with the most precise support
        long preciseFirst = first.support().stream().filter(x -> recordSoftClipsAtJunction(x.cachedRead(), first.junction())).count();
        long preciseSecond = second.support().stream().filter(x -> recordSoftClipsAtJunction(x.cachedRead(), second.junction())).count();

        return preciseFirst >= preciseSecond;
    }
}
