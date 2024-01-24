package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MERGE_READ_SUPPORT_OVERLAP;
import static com.hartwig.hmftools.esvee.common.AssemblyUtils.basesMatch;
import static com.hartwig.hmftools.esvee.common.AssemblyUtils.haveOverlappingReads;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import com.hartwig.hmftools.esvee.common.JunctionAssembly;

public class AssemblyDeduper
{
    public static void dedupJunctionAssemblies(final List<JunctionAssembly> assemblies)
    {
        dedupByAssemblyContainsAnother(assemblies);

    }

    private static void dedupByAssemblyContainsAnother(final List<JunctionAssembly> assemblies)
    {
        Collections.sort(assemblies, Collections.reverseOrder(Comparator.comparingInt(x -> x.length())));

        int i = 0;
        while(i < assemblies.size())
        {
            JunctionAssembly first = assemblies.get(i);

            int j = i + 1;
            while(j < assemblies.size())
            {
                JunctionAssembly second = assemblies.get(j);

                if(assemblyContainsAnother(first, second))
                {
                    assemblies.remove(j);
                    first.checkAddReadSupport(second);

                    continue;
                }

                ++j;
            }

            ++i;
        }
    }

    private static boolean assemblyContainsAnother(final JunctionAssembly first, final JunctionAssembly second)
    {
        int rangeStart;
        int rangeEnd;

        if(first.initialJunction().isForward())
        {
            rangeStart = first.initialJunction().Position;
            rangeEnd = min(first.maxAlignedPosition(), second.maxAlignedPosition());
        }
        else
        {
            rangeStart = max(first.minAlignedPosition(), second.minAlignedPosition());
            rangeEnd = first.initialJunction().Position;
        }

        int firstOffset = rangeStart - first.minAlignedPosition();;
        int secondOffset = rangeStart - second.minAlignedPosition();
        int rangeLength = rangeEnd - rangeStart;

        for(int i = 0; i <= rangeLength; ++i)
        {
            if(!basesMatch(
                    first.bases()[i + firstOffset], second.bases()[i + secondOffset],
                    first.baseQuals()[i + firstOffset], second.baseQuals()[i + secondOffset], LOW_BASE_QUAL_THRESHOLD))
            {
                return false;
            }
        }

        return true;
    }

}
