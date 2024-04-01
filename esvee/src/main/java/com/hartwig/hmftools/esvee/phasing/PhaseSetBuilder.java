package com.hartwig.hmftools.esvee.phasing;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PROXIMATE_DUP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.isLocalAssemblyCandidate;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.isSupplementaryOnly;
import static com.hartwig.hmftools.esvee.assembly.RemoteRegionAssembler.assemblyOverlapsRemoteRegion;
import static com.hartwig.hmftools.esvee.types.AssemblyOutcome.DUP_BRANCHED;
import static com.hartwig.hmftools.esvee.types.AssemblyOutcome.LINKED;
import static com.hartwig.hmftools.esvee.assembly.RefBaseExtender.extendRefBases;
import static com.hartwig.hmftools.esvee.assembly.AssemblyLinker.tryAssemblyFacing;
import static com.hartwig.hmftools.esvee.types.AssemblyOutcome.REMOTE_REF;
import static com.hartwig.hmftools.esvee.types.AssemblyOutcome.SHORT_INDEL;
import static com.hartwig.hmftools.esvee.types.AssemblySupport.findMatchingFragmentSupport;
import static com.hartwig.hmftools.esvee.types.AssemblySupport.hasMatchingFragment;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.assembliesShareReads;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.assembly.AssemblyLinker;
import com.hartwig.hmftools.esvee.assembly.IndelBuilder;
import com.hartwig.hmftools.esvee.assembly.RemoteRegionAssembler;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.types.AssemblyLink;
import com.hartwig.hmftools.esvee.types.AssemblySupport;
import com.hartwig.hmftools.esvee.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.types.LinkType;
import com.hartwig.hmftools.esvee.types.PhaseGroup;
import com.hartwig.hmftools.esvee.types.PhaseSet;
import com.hartwig.hmftools.esvee.types.RemoteRegion;
import com.hartwig.hmftools.esvee.types.SupportType;

public class PhaseSetBuilder
{
    private final PhaseGroup mPhaseGroup;
    private final RefGenomeInterface mRefGenome;
    private final RemoteRegionAssembler mRemoteRegionAssembler;

    // references from phase group
    private final List<JunctionAssembly> mAssemblies;
    private final List<PhaseSet> mPhaseSets; // final proposed phase sets

    // working cache only
    private final List<AssemblyLink> mSplitLinks;
    private final List<AssemblyLink> mFacingLinks;
    private final List<AssemblyLink> mSecondarySplitLinks;

    public PhaseSetBuilder(
            final RefGenomeInterface refGenome, final RemoteRegionAssembler remoteRegionAssembler, final PhaseGroup phaseGroup)
    {
        mRefGenome = refGenome;
        mPhaseGroup = phaseGroup;
        mRemoteRegionAssembler = remoteRegionAssembler;
        mPhaseSets = mPhaseGroup.phaseSets();
        mAssemblies = mPhaseGroup.assemblies();
        mSecondarySplitLinks = mPhaseGroup.secondaryLinks();

        mSplitLinks = Lists.newArrayList();
        mFacingLinks = Lists.newArrayList();
    }

    public void buildPhaseSets()
    {
        markShortIndelAssemblies();

        if(mAssemblies.size() == 2)
        {
            handleAssemblyPair();
            return;
        }

        formSplitLinks();

        extendRemoteAssemblies();

        formFacingLinks();

        // extendJunctions();

        formPhaseSets();

        cleanupBranchedAssemblies();
    }

    private void handleAssemblyPair()
    {
        // simpler routine without prioritising pairs, facing links or branching
        JunctionAssembly assembly1 = mAssemblies.get(0);
        JunctionAssembly assembly2 = mAssemblies.get(1);

        if(!isLocalAssemblyCandidate(assembly1, assembly2) && (assembly1.outcome() == SHORT_INDEL) || assembly2.outcome() == SHORT_INDEL)
            return;

        // if no link was made, then may need to revert to logic for finding discordant pair assemblies etc
        // likewise may depend for solo-assemblies on how disc-only assemblies are used
        AssemblyLink assemblyLink = checkSplitLink(assembly1, assembly2);

        if(assemblyLink != null)
        {
            buildSplitLink(assembly1, assembly2, false);
            mPhaseSets.add(new PhaseSet(assemblyLink));
        }
    }

    private void markShortIndelAssemblies()
    {
        // look for evidence that a local assembly with soft-clipped support is explained by a local short INDEL (length 10-31)
        for(JunctionAssembly assembly : mAssemblies)
        {
            if(assembly.indel())
                continue;

            boolean hasSupportingIndel = assembly.support().stream()
                    .anyMatch(x -> x.type() == SupportType.JUNCTION && IndelBuilder.convertedIndelCrossesJunction(assembly, x.read()));

            if(hasSupportingIndel)
            {
                assembly.setOutcome(SHORT_INDEL);
            }
        }
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

                boolean isLocalLink = isLocalAssemblyCandidate(assembly1, assembly2);

                if(!isLocalLink && (assembly1.outcome() == SHORT_INDEL || assembly2.outcome() == SHORT_INDEL))
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
                boolean allowBranching = !(assemblyLink.svType() == DUP && assemblyLink.length() < PROXIMATE_DUP_LENGTH);
                buildSplitLink(assembly1, assembly2, allowBranching);
                continue;
            }

            // form secondaries for any assemblies which aren't duplicates of another or supplementary only
            if(assembly1.outcome().isDuplicate() || assembly2.outcome().isDuplicate())
                continue;

            if(isSupplementaryOnly(assembly1) || isSupplementaryOnly(assembly2))
                continue;

            AssemblyLink assemblyLink = AssemblyLinker.tryAssemblyOverlap(assembly1, assembly2, true);

            if(assemblyLink != null)
                mSecondarySplitLinks.add(assemblyLink);
        }
    }

    private AssemblyLink checkSplitLink(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        if(assembly1.junction() == assembly2.junction()) // ignore duplicates
            return null;

        // handle local INDELs here since the following logic currently applies to them
        AssemblyLink assemblyLink = AssemblyLinker.tryAssemblyIndel(assembly1, assembly2);

        if(assemblyLink != null)
            return assemblyLink;

        return AssemblyLinker.tryAssemblyOverlap(assembly1, assembly2);
    }

    private void buildSplitLink(final JunctionAssembly assembly1, final JunctionAssembly assembly2, boolean allowBranching)
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
        extendRefBases(assembly1, matchedCandidates1, mRefGenome, allowBranching);
        extendRefBases(assembly2, matchedCandidates2, mRefGenome, allowBranching);
    }

    private static void checkMatchingCandidateSupport(
            final JunctionAssembly otherAssembly,
            final List<AssemblySupport> candidateSupport, final List<AssemblySupport> otherCandidateSupport,
            final List<AssemblySupport> matchedCandidates, final List<AssemblySupport> otherMatchedCandidates)
    {
        // consider each candidate support read to see if it has a matching read in the other assembly's candidates or junction reads
        for(AssemblySupport support : candidateSupport)
        {
            // add any junction mate reads local to the assembly
            if(support.type() == SupportType.JUNCTION_MATE)
            {
                matchedCandidates.add(support);
                continue;
            }

            // now check for discordant reads with matching support in the other assembly
            List<AssemblySupport> matchedSupport = findMatchingFragmentSupport(otherAssembly.support(), support.read());

            List<AssemblySupport> matchedCandidateSupport = findMatchingFragmentSupport(otherCandidateSupport, support.read());

            // remove from other's candidates to avoid checking again
            matchedCandidateSupport.forEach(x -> otherCandidateSupport.remove(x));
            otherMatchedCandidates.addAll(matchedCandidateSupport);

            matchedSupport.addAll(matchedCandidateSupport);

            if(!matchedSupport.isEmpty())
                matchedCandidates.add(support);

            for(AssemblySupport matched : matchedSupport)
            {
                support.read().makeReadLinks(matched.read());
            }
        }
    }

    private static void linkExistingSupport(final JunctionAssembly first, final JunctionAssembly second)
    {
        for(AssemblySupport support : first.support())
        {
            List<AssemblySupport> matchedSupport = findMatchingFragmentSupport(second.support(), support.read());
            matchedSupport.forEach(x -> support.read().makeReadLinks(x.read()));
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

    private void extendRemoteAssemblies()
    {
        // limit this to a subset of unliked assemblies:
        // sufficient evidence and quality, and with remote junction mates
        List<JunctionAssembly> unlinkedAssemblies = mAssemblies.stream()
                .filter(x -> RemoteRegionAssembler.isExtensionCandidateAssembly(x))
                .filter(x -> mSplitLinks.stream().noneMatch(y -> y.hasAssembly(x))).collect(Collectors.toList());

        boolean foundRemoteLink = false;

        for(JunctionAssembly assembly : unlinkedAssemblies)
        {
            // collect remote regions which aren't only supplementaries nor which overlap another phase assembly
            List<RemoteRegion> remoteRegions = assembly.remoteRegions().stream()
                    .filter(x -> !x.isSuppOnlyRegion())
                    .filter(x -> mAssemblies.stream().filter(y -> y != assembly).noneMatch(y -> assemblyOverlapsRemoteRegion(y, x)))
                    .collect(Collectors.toList());

            if(remoteRegions.isEmpty())
                continue;

            // evaluate by remote regions with most linked reads
            Collections.sort(remoteRegions, Comparator.comparingInt(x -> -x.nonSuppReadCount()));

            for(RemoteRegion remoteRegion : remoteRegions)
            {
                // could take candidate discordant reads that map to the same region
                List<Read> localReads = assembly.support().stream()
                        .filter(x -> remoteRegion.readIds().contains(x.read().id()))
                        .map(x -> x.read())
                        .collect(Collectors.toList());

                assembly.candidateSupport().stream()
                        .filter(x -> x.type() == SupportType.CANDIDATE_DISCORDANT)
                        .filter(x -> remoteRegion.readIds().contains(x.read().id()))
                        .forEach(x -> localReads.add(x.read()));

                AssemblyLink assemblyLink = mRemoteRegionAssembler.tryRemoteAssemblyLink(assembly, remoteRegion, localReads);

                if(assemblyLink == null)
                    continue;

                JunctionAssembly remoteAssembly = assemblyLink.otherAssembly(assembly);

                remoteAssembly.setOutcome(REMOTE_REF);

                mPhaseGroup.addDerivedAssembly(remoteAssembly);

                if(!foundRemoteLink)
                {
                    foundRemoteLink = true;
                    mSplitLinks.add(assemblyLink);

                    buildSplitLink(assembly, remoteAssembly, true);
                }
                else
                {
                    mSecondarySplitLinks.add(assemblyLink);
                }
            }
        }

        /*
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
        */
    }

    private void formFacingLinks()
    {
        if(mAssemblies.size() <= 2)
            return;

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
                        linkExistingSupport(assembly, facingAssembly);
                        facingAssemblies.add(facingAssembly);
                        facingAssemblies.add(assembly);
                    }
                }
            }
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

        for(PhaseSet phaseSet : mPhaseSets)
        {
            phaseSet.assemblies().forEach(x -> x.setOutcome(LINKED));
        }
    }

    private AssemblyLink findLinkedAssembly(final JunctionAssembly assembly, boolean findSplit)
    {
        // find a link using one assembly of a particular type, then remove it from future consideration
        List<AssemblyLink> searchLinks = findSplit ? mSplitLinks : mFacingLinks;

        int index = 0;
        while(index < searchLinks.size())
        {
            AssemblyLink link = searchLinks.get(index);

            if(link.hasAssembly(assembly))
            {
                searchLinks.remove(index);

                if(!findSplit)
                {
                    // remove any other facing links which use this assembly
                    JunctionAssembly otherAssembly = link.otherAssembly(assembly);

                    int otherIndex = 0;
                    while(otherIndex < mFacingLinks.size())
                    {
                        AssemblyLink otherLink = searchLinks.get(otherIndex);
                        if(otherLink.hasAssembly(assembly) || otherLink.hasAssembly(otherAssembly))
                            searchLinks.remove(otherLink);
                        else
                            ++otherIndex;
                    }
                }

                return link;
            }

            ++index;
        }

        return null;
    }
    private void cleanupBranchedAssemblies()
    {
        // finally remove any branched assemblies which did not form a facing link
        List<JunctionAssembly> branchedAssemblies = mAssemblies.stream()
                .filter(x -> x.outcome() == DUP_BRANCHED)
                .collect(Collectors.toList());

        for(JunctionAssembly branchedAssembly : branchedAssemblies)
        {
            boolean phaseLinked = false;

            for(PhaseSet phaseSet : mPhaseSets)
            {
                if(phaseSet.assemblyLinks().stream().filter(x -> x.type() == LinkType.FACING).anyMatch(x -> x.hasAssembly(branchedAssembly)))
                {
                    phaseLinked = true;
                    break;
                }
            }

            if(!phaseLinked)
            {
                mPhaseGroup.assemblies().remove(branchedAssembly);
                mPhaseGroup.derivedAssemblies().remove(branchedAssembly);
            }
        }
    }
}
