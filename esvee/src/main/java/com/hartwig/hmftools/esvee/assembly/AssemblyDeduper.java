package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_DEDUP_HIGH_SUPPORT_RATIO;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_DEDUP_JITTER_MAX_DIST;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PROXIMATE_JUNCTION_DISTANCE;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.indelAssembliesMatch;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.recordSoftClipsAtJunction;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;
import com.hartwig.hmftools.esvee.common.IndelCoords;

public class AssemblyDeduper
{
    public static void dedupProximateAssemblies(
            final List<JunctionAssembly> existingAssemblies, final List<JunctionAssembly> newAssemblies,
            final List<JunctionAssembly> dedupedIndels)
    {
        if(newAssemblies.isEmpty() || existingAssemblies.isEmpty())
            return;

        // start with the most recent previous assemblies since they are added in order
        int index = existingAssemblies.size() - 1;
        int minPosition = newAssemblies.get(0).junction().Position - PROXIMATE_JUNCTION_DISTANCE;

        List<JunctionAssembly> removeExisting = Lists.newArrayList();

        while(index >= 0 && !newAssemblies.isEmpty())
        {
            JunctionAssembly assembly = existingAssemblies.get(index);

            if(assembly.junction().Position < minPosition)
                break;

            int newIndex = 0;

            while(newIndex < newAssemblies.size())
            {
                JunctionAssembly newAssembly = newAssemblies.get(newIndex);

                if(newAssembly.junction() == assembly.junction() || newAssembly.junction().Orient != assembly.junction().Orient)
                {
                    ++newIndex;
                    continue;
                }

                if(newAssembly.indel() && assembly.indel())
                {
                    if(checkDedupIndelPair(newAssembly, assembly, newIndex, newAssemblies, existingAssemblies, removeExisting, dedupedIndels))
                        ++newIndex;
                }
                else
                {
                    if(checkDedupNonIndel(newAssembly, assembly, newIndex, newAssemblies, removeExisting))
                        ++newIndex;
                }
            }

            --index;
        }

        removeExisting.forEach(x -> existingAssemblies.remove(x));
    }

    private static boolean checkDedupNonIndel(
            final JunctionAssembly newAssembly, final JunctionAssembly assembly, final int newIndex,
            final List<JunctionAssembly> newAssemblies, final List<JunctionAssembly> removeExisting)
    {
        // removes true if keeps the new assembly
        if(!canDedupAssemblies(assembly, newAssembly))
            return true;

        if(selectFirstAssembly(assembly, newAssembly))
        {
            assembly.addMergedAssembly();
            newAssemblies.remove(newIndex);
            return false;
        }
        else
        {
            // keep this assembly and mark the existing one for removal
            newAssembly.addMergedAssembly();
            removeExisting.add(assembly);
            return true;
        }
    }

    private static boolean checkDedupIndelPair(
            final JunctionAssembly newAssembly, final JunctionAssembly assembly, final int newIndex,
            final List<JunctionAssembly> newAssemblies, final List<JunctionAssembly> existingAssemblies,
            final List<JunctionAssembly> removeExisting, final List<JunctionAssembly> dedupedIndels)
    {
        // removes true if keeps the new assembly
        JunctionAssembly matchedIndel = newAssembly.indel() ? dedupedIndels.stream()
                .filter(x -> indelAssembliesMatch(x, newAssembly)).findFirst().orElse(null) : null;

        boolean dedupNewAssembly = false;

        if(matchedIndel == null)
        {
            if(!canDedupAssemblies(assembly, newAssembly))
                return true;

            // take the assembly with the most read support
            dedupNewAssembly = selectFirstAssembly(assembly, newAssembly);
        }
        else
        {
            dedupNewAssembly = true;
            dedupedIndels.remove(matchedIndel);
        }

        if(dedupNewAssembly)
        {
            assembly.addMergedAssembly();
            newAssemblies.remove(newIndex);

            if(matchedIndel == null)
            {
                // remove the paired other assembly, as long the retained assembly has a matching pair
                JunctionAssembly matchedExistingIndel = existingAssemblies.stream()
                        .filter(x -> indelAssembliesMatch(x, newAssembly)).findFirst().orElse(null);

                if(matchedExistingIndel != null)
                {
                    if(existingAssemblies.stream().anyMatch(x -> indelAssembliesMatch(x, assembly)))
                        removeExisting.add(matchedExistingIndel);
                }
                else
                {
                    dedupedIndels.add(newAssembly);
                }
            }

            return false;
        }
        else
        {
            // keep this assembly and mark the existing one for removal
            newAssembly.addMergedAssembly();
            removeExisting.add(assembly);

            JunctionAssembly matchedExistingIndel = existingAssemblies.stream()
                    .filter(x -> indelAssembliesMatch(x, assembly)).findFirst().orElse(null);

            if(matchedExistingIndel != null)
            {
                if(existingAssemblies.stream().anyMatch(x -> indelAssembliesMatch(x, newAssembly)))
                    removeExisting.add(matchedExistingIndel);
            }

            return true;
        }
    }

    private static boolean canDedupAssemblies(final JunctionAssembly first, final JunctionAssembly second)
    {
        if(dedupOnSupportRatio(first, second))
            return true;

        if(!SequenceCompare.matchedAssemblySequences(first, second))
            return false;

        if(first.indel() == second.indel()) // remaining checks are for indel vs non-indel
            return true;

        // despite a sequence match, if an indel read supports a soft-clip junction without offering inferred support, do not dedup
        int junctionDistance = abs(first.junction().Position - second.junction().Position);

        if(junctionDistance <= ASSEMBLY_DEDUP_JITTER_MAX_DIST) // must be sufficiently far apart
            return true;

        IndelCoords firstIndelCoords = first.indelCoords();
        IndelCoords secondIndelCoords = second.indelCoords();

        if(secondIndelCoords != null)
        {
            for(SupportRead firstRead : first.support())
            {
                if(firstRead.type() != SupportType.JUNCTION)
                    continue;

                if(firstRead.indelCoords() != null && firstRead.indelCoords().matches(secondIndelCoords))
                    return false;
            }
        }

        if(firstIndelCoords != null)
        {
            for(SupportRead secondRead : second.support())
            {
                if(secondRead.type() != SupportType.JUNCTION)
                    continue;

                if(secondRead.indelCoords() != null && secondRead.indelCoords().matches(firstIndelCoords))
                    return false;
            }
        }

        return true;
    }

    private static boolean dedupOnSupportRatio(final JunctionAssembly first, final JunctionAssembly second)
    {
        // dedup an assembly if its support is much less than the other assembly
        int junctionDistance = abs(first.junction().Position - second.junction().Position);

        if(junctionDistance > ASSEMBLY_DEDUP_JITTER_MAX_DIST)
            return false;

        return assemblyExceedsSupportRatio(first, second) || assemblyExceedsSupportRatio(second, first);
    }

    private static boolean assemblyExceedsSupportRatio(final JunctionAssembly first, final JunctionAssembly second)
    {
        return first.supportCount() > second.supportCount() * ASSEMBLY_DEDUP_HIGH_SUPPORT_RATIO;
    }

    private static boolean selectFirstAssembly(final JunctionAssembly first, final JunctionAssembly second)
    {
        if(assemblyExceedsSupportRatio(first, second))
            return true;

        if(assemblyExceedsSupportRatio(second, first))
            return false;

        if(first.indel() != second.indel())
            return first.indel();

        // always favour assemblies with germline over those without, so as not to miss germline support
        boolean firstHasGermlineSupport = first.support().stream().anyMatch(x -> x.isReference());
        boolean secondHasGermlineSupport = second.support().stream().anyMatch(x -> x.isReference());

        if(firstHasGermlineSupport != secondHasGermlineSupport)
            return firstHasGermlineSupport;

        if(first.discordantOnly() != second.discordantOnly())
            return !first.discordantOnly();

        if(first.supportCount() != second.supportCount())
            return first.supportCount() > second.supportCount();

        // take the one with the most precise support
        long preciseFirst = first.support().stream().filter(x -> recordSoftClipsAtJunction(x.cachedRead(), first.junction())).count();
        long preciseSecond = second.support().stream().filter(x -> recordSoftClipsAtJunction(x.cachedRead(), second.junction())).count();

        return preciseFirst >= preciseSecond;
    }
}
