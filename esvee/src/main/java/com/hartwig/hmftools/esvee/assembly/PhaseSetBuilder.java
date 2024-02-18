package com.hartwig.hmftools.esvee.assembly;

import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_LENGTH;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyExtender.extendRefBases;
import static com.hartwig.hmftools.esvee.common.AssemblySupport.findMatchingFragmentSupport;
import static com.hartwig.hmftools.esvee.common.AssemblySupport.hasMatchingFragment;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.common.AssemblyLink;
import com.hartwig.hmftools.esvee.common.AssemblySupport;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.PhaseGroup;
import com.hartwig.hmftools.esvee.common.PhaseSet;
import com.hartwig.hmftools.esvee.common.RefBaseAssembly;
import com.hartwig.hmftools.esvee.common.RefSideSoftClip;

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
        formSplitLinks();

        for(JunctionAssembly assembly : mPhaseGroup.assemblies())
        {
            for(JunctionAssembly branchedAssembly : assembly.branchedAssemblies())
            {
                mPhaseGroup.addAssembly(branchedAssembly);

                // replicate the links? or model how?
            }
        }

        formFacingLinks();
    }

    private void formFacingLinks()
    {
        List<JunctionAssembly> assemblies = mPhaseGroup.assemblies();

        if(assemblies.size() <= 2)
            return;

        // support will have changed now, so reassess candidates for facing links

    }

    private void formSplitLinks()
    {
        List<JunctionAssembly> assemblies = mPhaseGroup.assemblies();

        if(assemblies.size() == 2)
        {
            // no need to re-test shared read support since they must be to be in a group
            AssemblyLink assemblyLink = checkSplitLink(assemblies.get(0), assemblies.get(1));

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

        // build any split links and only allow an assembly to be used once
        Set<JunctionAssembly> linkedAssemblies = Sets.newHashSet();

        while(!assemblySupportPairs.isEmpty())
        {
            SharedAssemblySupport sharedReadPair = assemblySupportPairs.get(0);
            assemblySupportPairs.remove(0);

            if(linkedAssemblies.contains(sharedReadPair.Assembly1) || linkedAssemblies.contains(sharedReadPair.Assembly2))
                continue;

            // test if a link can be made
            AssemblyLink splitLink = checkSplitLink(sharedReadPair.Assembly1, sharedReadPair.Assembly2);

            if(splitLink == null)
                continue;

            linkedAssemblies.add(sharedReadPair.Assembly1);
            linkedAssemblies.add(sharedReadPair.Assembly2);
        }

        /*
        while(!assemblySupportPairs.isEmpty())
        {
            SharedAssemblySupport sharedReadPair = assemblySupportPairs.get(0);
            assemblySupportPairs.remove(0);

            // skip if already matches a link
            if(assemblyLinkExists(sharedReadPair.Assembly1, sharedReadPair.Assembly2))
                continue;

            // test if a link can be made
            AssemblyLink sharedReadLink = checkSplitLink(sharedReadPair.Assembly1, sharedReadPair.Assembly2);

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
        */
    }

    private AssemblyLink checkSplitLink(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        AssemblyLink assemblyLink = mAssemblyLinker.tryAssemblyOverlap(assembly1, assembly2);

        if(assemblyLink == null)
            return null;

        // FIXME: make links between junction reads

        // look for shared reads between the assemblies, and factor in discordant reads which were only considered candidates until now
        List<AssemblySupport> matchedCandidates1 = Lists.newArrayList();
        List<AssemblySupport> matchedCandidates2 = Lists.newArrayList();

        // list is copied since matched candidates will be removed from repeated matching
        List<AssemblySupport> candidateSupport2 = Lists.newArrayList(assembly2.candidateSupport());

        checkMatchingCandidateSupport(assembly2, assembly1.candidateSupport(), candidateSupport2, matchedCandidates1, matchedCandidates2);
        checkMatchingCandidateSupport(assembly1, candidateSupport2, Collections.emptyList(), matchedCandidates2, matchedCandidates1);

        // build out ref-base assembly support from this non-junction reads
        extendRefBases(assembly1, matchedCandidates1);
        extendRefBases(assembly1, matchedCandidates2);

        return assemblyLink;
    }

    private static void checkMatchingCandidateSupport(
            final JunctionAssembly otherAssembly,
            final List<AssemblySupport> candidateSupport, final List<AssemblySupport> otherCandidateSupport,
            final List<AssemblySupport> matchedCandidates, final List<AssemblySupport> otherMatchedCandidates)
    {
        // TODO: at this point reads could be linked to each other - any reason not to do this?
        //  could make fragment count determination more efficient
        for(AssemblySupport support : candidateSupport)
        {
            AssemblySupport matchedSupport = findMatchingFragmentSupport(otherAssembly.support(), support.read());

            if(matchedSupport == null)
            {
                matchedSupport = findMatchingFragmentSupport(otherCandidateSupport, support.read());

                if(matchedSupport != null)
                {
                    // remove from other's candidates to avoid checking again
                    otherMatchedCandidates.add(matchedSupport);
                    otherCandidateSupport.remove(matchedSupport);
                }
            }

            if(matchedSupport != null)
            {
                matchedCandidates.add(support);
                support.read().makeReadLinks(matchedSupport.read());
            }
        }
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

                    AssemblyLink phasedLink = checkSplitLink(newAssembly, existingAssembly);

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
