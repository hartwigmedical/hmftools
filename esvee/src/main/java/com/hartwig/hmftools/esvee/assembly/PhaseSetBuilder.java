package com.hartwig.hmftools.esvee.assembly;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.esvee.assembly.AssemblyOutcome.DUP_SPLIT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyOutcome.LINKED;
import static com.hartwig.hmftools.esvee.assembly.PhaseGroupBuilder.isLocalAssemblyCandidate;
import static com.hartwig.hmftools.esvee.assembly.RefBaseExtender.extendRefBases;
import static com.hartwig.hmftools.esvee.assembly.AssemblyLinker.tryAssemblyFacing;
import static com.hartwig.hmftools.esvee.common.AssemblySupport.findMatchingFragmentSupport;
import static com.hartwig.hmftools.esvee.common.AssemblySupport.hasMatchingFragment;
import static com.hartwig.hmftools.esvee.common.AssemblyUtils.assembliesShareReads;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.common.AssemblyLink;
import com.hartwig.hmftools.esvee.common.AssemblySupport;
import com.hartwig.hmftools.esvee.common.DiscordantGroup;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.PhaseGroup;
import com.hartwig.hmftools.esvee.common.PhaseSet;
import com.hartwig.hmftools.esvee.common.SupportType;

public class PhaseSetBuilder
{
    private final PhaseGroup mPhaseGroup;
    private final RefGenomeInterface mRefGenome;

    private final AssemblyLinker mAssemblyLinker;

    // references from phase group
    private final List<JunctionAssembly> mAssemblies;
    private final List<PhaseSet> mPhaseSets; // final proposed phase sets

    // working cache only
    private final List<AssemblyLink> mSplitLinks;
    private final List<AssemblyLink> mFacingLinks;
    private final List<AssemblyLink> mSecondarySplitLinks;

    public PhaseSetBuilder(final RefGenomeInterface refGenome, final PhaseGroup phaseGroup)
    {
        mRefGenome = refGenome;
        mPhaseGroup = phaseGroup;
        mAssemblyLinker = new AssemblyLinker();
        mPhaseSets = mPhaseGroup.phaseSets();
        mAssemblies = mPhaseGroup.assemblies();
        mSecondarySplitLinks = mPhaseGroup.secondaryLinks();

        mSplitLinks = Lists.newArrayList();
        mFacingLinks = Lists.newArrayList();
    }

    public void buildPhaseSets()
    {
        if(mAssemblies.size() == 2)
        {
            // if no link was made, then may need to revert to logic for finding discordant pair assemblies etc
            // likewise may depend for solo-assemblies how disc-only assemblies are used
            AssemblyLink assemblyLink = checkSplitLink(mAssemblies.get(0), mAssemblies.get(1));

            if(assemblyLink != null)
            {
                buildSplitLink(mAssemblies.get(0), mAssemblies.get(1));
                mPhaseSets.add(new PhaseSet(assemblyLink));
            }

            return;
        }

        formSplitLinks();

        // add any branched assemblies to the phase group
        List<JunctionAssembly> branchedAssemblies = Lists.newArrayList();
        mPhaseGroup.assemblies().forEach(x -> branchedAssemblies.addAll(x.branchedAssemblies()));
        branchedAssemblies.forEach(x -> mPhaseGroup.addAssembly(x));

        formFacingLinks();

        // extendJunctions();

        formPhaseSets();
    }

    private void formSplitLinks()
    {
        // where there are more than 2 assemblies, start with the ones with the most support and overlapping junction reads
        List<SharedAssemblySupport> assemblySupportPairs = Lists.newArrayList();

        for(int i = 0; i < mAssemblies.size() - 1; ++i)
        {
            JunctionAssembly assembly1 = mAssemblies.get(i);

            for(int j = i + 1; j < mAssemblies.size(); ++j)
            {
                JunctionAssembly assembly2 = mAssemblies.get(j);

                if(assembly1.indel() != assembly2.indel())
                    continue;

                int sharedCount = 0;

                for(AssemblySupport support : assembly1.support())
                {
                    if(hasMatchingFragment(assembly2.support(), support.read()))
                        ++sharedCount;
                }

                // also check candidate reads
                for(AssemblySupport support : assembly1.candidateSupport())
                {
                    if(hasMatchingFragment(assembly2.support(), support.read()))
                        ++sharedCount;
                }

                boolean isLocalLink = isLocalAssemblyCandidate(assembly1, assembly2);

                if(sharedCount > 1 || isLocalLink)
                    assemblySupportPairs.add(new SharedAssemblySupport(assembly1, assembly2, sharedCount, isLocalLink));
            }
        }

        Collections.sort(assemblySupportPairs);

        // build any split links and only allow an assembly to be used once
        Set<JunctionAssembly> linkedAssemblies = Sets.newHashSet();

        while(!assemblySupportPairs.isEmpty())
        {
            SharedAssemblySupport sharedReadPair = assemblySupportPairs.remove(0);

            JunctionAssembly assembly1 = sharedReadPair.Assembly1;
            JunctionAssembly assembly2 = sharedReadPair.Assembly2;

            boolean alreadyLinked = linkedAssemblies.contains(assembly1) || linkedAssemblies.contains(assembly2);

            /*
            boolean hasLinkedDuplicate =
                    (assembly1.outcome().isDuplicate() && linkedAssemblies.stream().anyMatch(x -> x.junction() == assembly1.junction()))
                    || (assembly2.outcome().isDuplicate() && linkedAssemblies.stream().anyMatch(x -> x.junction() == assembly2.junction()));
            */

            // test if a link can be made
            if(!alreadyLinked)
            {
                AssemblyLink assemblyLink = checkSplitLink(assembly1, assembly2);

                if(assemblyLink == null)
                    continue;

                mSplitLinks.add(assemblyLink);
                linkedAssemblies.add(assembly1);
                linkedAssemblies.add(assembly2);
                buildSplitLink(assembly1, assembly2);
                continue;
            }

            // form secondaries for any assemblies which aren't duplicates of another
            if(assembly1.outcome().isDuplicate() || assembly2.outcome().isDuplicate())
                continue;

            AssemblyLink assemblyLink = mAssemblyLinker.tryAssemblyOverlap(assembly1, assembly2);

            if(assemblyLink != null)
                mSecondarySplitLinks.add(assemblyLink);
        }
    }

    private AssemblyLink checkSplitLink(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        if(assembly1.junction() == assembly2.junction()) // ignore duplicates
            return null;

        // handle local INDELs here since the following logic currently applies to them
        AssemblyLink assemblyLink = mAssemblyLinker.tryAssemblyIndel(assembly1, assembly2);

        if(assemblyLink != null)
            return assemblyLink;

        return mAssemblyLinker.tryAssemblyOverlap(assembly1, assembly2);
    }

    private void buildSplitLink(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        linkExistingSupport(assembly1, assembly2);

        assembly1.setOutcome(LINKED);
        assembly2.setOutcome(LINKED);

        // look for shared reads between the assemblies, and factor in discordant reads which were only considered candidates until now
        List<AssemblySupport> matchedCandidates1 = Lists.newArrayList();
        List<AssemblySupport> matchedCandidates2 = Lists.newArrayList();

        // list is copied since matched candidates will be removed from repeated matching
        List<AssemblySupport> candidateSupport2 = Lists.newArrayList(assembly2.candidateSupport());

        // find matching reads, and link reads to each other where possible
        checkMatchingCandidateSupport(assembly2, assembly1.candidateSupport(), candidateSupport2, matchedCandidates1, matchedCandidates2);
        checkMatchingCandidateSupport(assembly1, candidateSupport2, Collections.emptyList(), matchedCandidates2, matchedCandidates1);

        // build out ref-base assembly support from these non-junction reads - both matched discordant and junction mates
        extendRefBases(assembly1, matchedCandidates1, mRefGenome);
        extendRefBases(assembly2, matchedCandidates2, mRefGenome);
    }

    private static void checkMatchingCandidateSupport(
            final JunctionAssembly otherAssembly,
            final List<AssemblySupport> candidateSupport, final List<AssemblySupport> otherCandidateSupport,
            final List<AssemblySupport> matchedCandidates, final List<AssemblySupport> otherMatchedCandidates)
    {
        for(AssemblySupport support : candidateSupport)
        {
            // add any junction mate reads local to the assembly
            if(support.type() == SupportType.JUNCTION_MATE)
            {
                matchedCandidates.add(support);
                continue;
            }

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

    private static void linkExistingSupport(final JunctionAssembly first, final JunctionAssembly second)
    {
        for(AssemblySupport support : first.support())
        {
            AssemblySupport matchedSupport = findMatchingFragmentSupport(second.support(), support.read());

            if(matchedSupport != null)
            {
                support.read().makeReadLinks(matchedSupport.read());
            }
        }
    }

    private class SharedAssemblySupport implements Comparable<SharedAssemblySupport>
    {
        public final JunctionAssembly Assembly1;
        public final JunctionAssembly Assembly2;
        public final int SharedSupport;
        public final boolean IsLocalDelDup;

        public SharedAssemblySupport(
                final JunctionAssembly assembly1, final JunctionAssembly assembly2, final int sharedSupport, final boolean localDelDup)
        {
            Assembly1 = assembly1;
            Assembly2 = assembly2;
            SharedSupport = sharedSupport;
            IsLocalDelDup = localDelDup;
        }

        @Override
        public int compareTo(final SharedAssemblySupport other)
        {
            if(other.IsLocalDelDup != IsLocalDelDup)
                return IsLocalDelDup ? -1 : 1;

            if(SharedSupport != other.SharedSupport)
                return SharedSupport > other.SharedSupport ? -1 : 1;

            return 0;
        }

        public String toString()
        {
            return format("%s + %s shared(%d) %s",
                    Assembly1.junction().coords(), Assembly2.junction().coords(), SharedSupport, IsLocalDelDup ? "localDelDup" : "");
        }
    }

    private void formFacingLinks()
    {
        if(mAssemblies.size() <= 2)
            return;

        // support will have changed now, so reassess candidates for facing links

        // for each assembly in a split link, look for a facing link (whether linked or not)
        Set<JunctionAssembly> facingAssemblies = Sets.newHashSet();

        for(AssemblyLink splitLink : mSplitLinks)
        {
            for(JunctionAssembly assembly : mAssemblies)
            {
                if(splitLink.hasAssembly(assembly) || facingAssemblies.contains(assembly))
                    continue;

                for(int i = 0; i <= 1; ++i)
                {
                    JunctionAssembly facingAssembly = (i == 0) ? splitLink.first() : splitLink.second();
                    JunctionAssembly splitAssembly = (i == 0) ? splitLink.second() : splitLink.first();

                    AssemblyLink facingLink = tryAssemblyFacing(facingAssembly, assembly);

                    if(facingLink == null)
                        continue;

                    // compelling evidence is a read from the new assembly which overlaps with the linked junction's reads
                    if(assembliesShareReads(assembly, splitAssembly))
                    {
                        mFacingLinks.add(facingLink);
                        facingAssemblies.add(facingAssembly);
                        facingAssemblies.add(assembly);
                    }
                }
            }
        }
    }

    private void extendJunctions()
    {
        // attempt to extend any assembly in a split link with unmapped mates, unlinked assemblies and discordant read groups
        List<JunctionAssembly> linkedAssemblies = Lists.newArrayList();
        List<JunctionAssembly> unlinkedAssemblies = Lists.newArrayList();

        for(JunctionAssembly assembly : mAssemblies)
        {
            if(mSplitLinks.stream().anyMatch(x -> x.hasAssembly(assembly)))
                linkedAssemblies.add(assembly);
            else
                unlinkedAssemblies.add(assembly);
        }

        List<DiscordantGroup> discordantGroups = mPhaseGroup.discordantGroups();

        for(JunctionAssembly assembly : linkedAssemblies)
        {
            JunctionExtender junctionExtender = new JunctionExtender(assembly, unlinkedAssemblies, discordantGroups);
            junctionExtender.extendAssembly();
        }
    }

    private void formPhaseSets()
    {
        // use split and facing links to assign assemblies to phase sets
        while(!mSplitLinks.isEmpty())
        {
            AssemblyLink splitLink = mSplitLinks.remove(0);

            PhaseSet phaseSet = new PhaseSet(splitLink);
            mPhaseSets.add(phaseSet);

            // look for facing and then splits links for this phase set
            for(int se = SE_START; se <= SE_END; ++se)
            {
                // check start and then end links of this phase set
                JunctionAssembly linkingAssembly = (se == SE_START) ? splitLink.first() : splitLink.second();
                boolean findSplit = false;

                while(true)
                {
                    AssemblyLink newLink = findLinkedAssembly(linkingAssembly, findSplit);

                    if(newLink == null)
                        break;

                    if(se == SE_START)
                        phaseSet.addAssemblyStart(newLink);
                    else
                        phaseSet.addAssemblyEnd(newLink);

                    findSplit = !findSplit;
                    linkingAssembly = newLink.otherAssembly(linkingAssembly);
                }
            }
        }
    }

    private AssemblyLink findLinkedAssembly(final JunctionAssembly assembly, boolean findSplit)
    {
        List<AssemblyLink> searchLinks = findSplit ? mSplitLinks : mFacingLinks;

        int index = 0;
        while(index < searchLinks.size())
        {
            AssemblyLink link = searchLinks.get(index);

            if(link.hasAssembly(assembly))
            {
                searchLinks.remove(index);
                return link;
            }

            ++index;
        }

        return null;
    }
}
