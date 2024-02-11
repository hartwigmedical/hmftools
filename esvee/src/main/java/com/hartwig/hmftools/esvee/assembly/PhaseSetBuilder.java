package com.hartwig.hmftools.esvee.assembly;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.PhaseGroup;
import com.hartwig.hmftools.esvee.common.PhaseSet;

public class PhaseSetBuilder
{
    private final PhaseGroup mPhaseGroup;

    private final List<PhaseSet> mPhaseSets;

    public PhaseSetBuilder(final PhaseGroup phaseGroup)
    {
        mPhaseGroup = phaseGroup;
        mPhaseSets = Lists.newArrayList();
    }

    public void buildPhaseSets()
    {
        List<JunctionAssembly> assemblies = mPhaseGroup.assemblies();

        if(assemblies.size() == 2)
        {
            // AssemblyLinker.tryAssemblyOverlap(assemblies.get(0), assemblies.get(1));
            return;
        }

        // where there are more than 2 assemblies, start with the ones with the most support and overlapping junction reads

        Collections.sort(assemblies, Comparator.comparingInt(x -> -x.supportCount()));

        List<JunctionAssembly> remainingAssemblies = assemblies.stream()
                .sorted(Comparator.comparingInt(x -> -x.supportCount()))
                .collect(Collectors.toList());

        // while(remainingAssemblies.isEmpty())


    }

}
