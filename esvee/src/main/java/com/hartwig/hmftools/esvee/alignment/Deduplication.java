package com.hartwig.hmftools.esvee.alignment;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.assembly.filters.FilterType;

public final class Deduplication
{
    private final List<Breakend> mDuplicateBreakends;

    public Deduplication()
    {
        mDuplicateBreakends = Lists.newArrayList();
    }

    public void deduplicateBreakends(final List<AssemblyAlignment> assemblyAlignments)
    {
        List<AssemblyBreakend> assemblyBreakends = Lists.newArrayList();

        for(AssemblyAlignment assemblyAlignment : assemblyAlignments)
        {
            for(Breakend breakend : assemblyAlignment.breakends())
            {
                assemblyBreakends.add(new AssemblyBreakend(assemblyAlignment, breakend));
            }
        }

        Collections.sort(assemblyBreakends);

        for(int i = 0; i < assemblyBreakends.size() - 1; ++i)
        {
            AssemblyBreakend assemblyBreakend = assemblyBreakends.get(i);

            if(!assemblyBreakend.BreakendRef.passing())
                continue;

            for(int j = i + 1; j < assemblyBreakends.size(); ++j)
            {
                AssemblyBreakend nextAssemblyBreakend = assemblyBreakends.get(j);

                // FIXME: adjust search based on position-matching rules
                if(isDuplicate(assemblyBreakend.BreakendRef, nextAssemblyBreakend.BreakendRef))
                {
                    // keep breakend with most support
                    if(keepFirst(assemblyBreakend.BreakendRef, nextAssemblyBreakend.BreakendRef))
                        nextAssemblyBreakend.BreakendRef.addFilter(FilterType.DUPLICATE);
                    else
                        assemblyBreakend.BreakendRef.addFilter(FilterType.DUPLICATE);
                }
                else
                {
                    break;
                }
            }
        }
    }

    private static boolean isDuplicate(final Breakend first, final Breakend second)
    {
        // CHECK: take confidence internal into consideration?
        return first.Chromosome.equals(second) && first.Position == second.Position && first.Orientation == second.Orientation;
    }

    private static boolean keepFirst(final Breakend first, final Breakend second)
    {
        int firstSupport = first.sampleSupport().stream().mapToInt(x -> x.totalSupport()).sum();
        int secondSupport = second.sampleSupport().stream().mapToInt(x -> x.totalSupport()).sum();
        return firstSupport >= secondSupport;
    }

    private class AssemblyBreakend implements Comparable<AssemblyBreakend>
    {
        public final AssemblyAlignment AssemblyAlignmentRef;
        public final Breakend BreakendRef;

        public AssemblyBreakend(final AssemblyAlignment assemblyAlignment, final Breakend breakend)
        {
            AssemblyAlignmentRef = assemblyAlignment;
            BreakendRef = breakend;
        }

        @Override
        public int compareTo(final AssemblyBreakend other)
        {
            return BreakendRef.compareTo(other.BreakendRef);
        }
    }
}
