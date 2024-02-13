package com.hartwig.hmftools.esvee.assembly;

import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.common.AssemblySupport.hasMatchingFragment;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.common.AssemblyLink;
import com.hartwig.hmftools.esvee.common.AssemblySupport;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.PhaseGroup;
import com.hartwig.hmftools.esvee.common.PhaseSet;

public class PhaseSetBuilder
{
    private final AssemblyLinker mAssemblyLinker;
    private final PhaseGroup mPhaseGroup;

    private final List<PhaseSet> mPhaseSets;

    public PhaseSetBuilder(final PhaseGroup phaseGroup)
    {
        mAssemblyLinker = new AssemblyLinker();
        mPhaseGroup = phaseGroup;
        mPhaseSets = mPhaseGroup.phaseSets();
    }

    public void buildPhaseSets()
    {
        List<JunctionAssembly> assemblies = mPhaseGroup.assemblies();

        if(assemblies.size() == 2)
        {
            // no need to re-test shared read support since they must be to be in a group
            AssemblyLink assemblyLink = mAssemblyLinker.tryAssemblyLink(assemblies.get(0), assemblies.get(1));

            if(assemblyLink != null)
            {
                mPhaseSets.add(new PhaseSet(assemblyLink));
            }

            return;
        }

        // where there are more than 2 assemblies, start with the ones with the most support and overlapping junction reads
        List<SharedAssemblySupport> assemblySupportPairs = Lists.newArrayList();

        for(int i = 0; i < assemblies.size() - 1; ++i)
        {
            JunctionAssembly assembly1 = assemblies.get(i);

            for(int j = i + 1; j < assemblies.size(); ++j)
            {
                JunctionAssembly assembly2 = assemblies.get(j);

                int sharedCount = 0;

                for(AssemblySupport support : assembly1.support())
                {
                    if(hasMatchingFragment(assembly2.support(), support.read()))
                        ++sharedCount;
                }

                if(sharedCount > 1)
                    assemblySupportPairs.add(new SharedAssemblySupport(assembly1, assembly2, sharedCount));
            }
        }

        Collections.sort(assemblySupportPairs, Comparator.comparingInt(x -> -x.SharedSupport));

        while(!assemblySupportPairs.isEmpty())
        {
            SharedAssemblySupport sharedReadPair = assemblySupportPairs.get(0);
            assemblySupportPairs.remove(0);

            // skip if already matches a link
            if(assemblyLinkExists(sharedReadPair.Assembly1, sharedReadPair.Assembly2))
                continue;

            // test if a link can be made
            AssemblyLink sharedReadLink = mAssemblyLinker.tryAssemblyLink(sharedReadPair.Assembly1, sharedReadPair.Assembly2);

            boolean addedToPhased = false;

            if(!mPhaseSets.isEmpty())
            {
                for(PhaseSet phaseSet : mPhaseSets)
                {
                    if(sharedReadLink != null
                    && (phaseSet.hasAssembly(sharedReadPair.Assembly1) || phaseSet.hasAssembly(sharedReadPair.Assembly2)))
                    {
                        phaseSet.addAssemblyEnd(sharedReadLink);
                        addedToPhased = true;
                        break;
                    }

                    AssemblyLink phasedLink = checkCanLinkWithPhaseSet(phaseSet, sharedReadPair);

                    if(phasedLink != null)
                    {
                        phaseSet.addAssemblyEnd(phasedLink);

                        if(sharedReadLink != null && !sharedReadLink.matches(phasedLink))
                            phaseSet.addAssemblyEnd(sharedReadLink);

                        addedToPhased = true;

                        // remove any entries of this kind from the shared-read-pair list to save reprocessing
                        for(int i = 0; i < assemblySupportPairs.size(); ++i)
                        {
                            SharedAssemblySupport next = assemblySupportPairs.get(i);
                            if(phasedLink.matches(next.Assembly1, next.Assembly2))
                            {
                                assemblySupportPairs.remove(i);
                                break;
                            }
                        }

                        break;
                    }
                }
            }

            if(!addedToPhased && sharedReadLink != null)
            {
                PhaseSet phaseSet = new PhaseSet(sharedReadLink);
                mPhaseSets.add(phaseSet);
            }
        }

        // TODO - merge any phase sets which can be and which share links? or no need?
    }

    private AssemblyLink checkCanLinkWithPhaseSet(final PhaseSet phaseSet, final SharedAssemblySupport sharedAssemblySupport)
    {
        // have already checked that no phase sets have this precise link, so now check if either of the new link's
        // assemblies can link with an existing assembly in this phase set
        for(AssemblyLink assemblyLink : phaseSet.assemblyLinks())
        {
            for(int n = 0; n <= 1; ++n)
            {
                JunctionAssembly newAssembly = (n == 0) ? sharedAssemblySupport.Assembly1 : sharedAssemblySupport.Assembly2;
                JunctionAssembly otherNewAssembly = (n == 0) ? sharedAssemblySupport.Assembly2 : sharedAssemblySupport.Assembly1;

                for(int e = 0; e <= 1; ++e)
                {
                    JunctionAssembly existingAssembly = (e == 0) ? assemblyLink.first() : assemblyLink.second();

                    if(existingAssembly == otherNewAssembly)
                        continue;

                    if(phaseSet.hasMatchingAssembly(newAssembly, existingAssembly))
                        return null;

                    AssemblyLink phasedLink = mAssemblyLinker.tryAssemblyLink(newAssembly, existingAssembly);

                    if(phasedLink != null)
                        return phasedLink;
                }
            }
        }

        return null;
    }

    private boolean assemblyLinkExists(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        return mPhaseSets.stream().anyMatch(x -> x.hasMatchingAssembly(assembly1, assembly2));
    }

    private class SharedAssemblySupport
    {
        public final JunctionAssembly Assembly1;
        public final JunctionAssembly Assembly2;
        public final int SharedSupport;

        public SharedAssemblySupport(final JunctionAssembly assembly1, final JunctionAssembly assembly2, final int sharedSupport)
        {
            Assembly1 = assembly1;
            Assembly2 = assembly2;
            SharedSupport = sharedSupport;
        }

        public String toString()
        {
            return format("%s + %s shared(%d)", Assembly1.junction().coords(), Assembly2.junction().coords(), SharedSupport);
        }
    }

}
