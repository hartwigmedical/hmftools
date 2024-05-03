package com.hartwig.hmftools.esvee.assembly.phase;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.esvee.AssemblyConstants.LOCAL_ASSEMBLY_MATCH_DISTANCE;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PROXIMATE_DUP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.isLocalAssemblyCandidate;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.isSupplementaryOnly;
import static com.hartwig.hmftools.esvee.assembly.RemoteRegionAssembler.assemblyOverlapsRemoteRegion;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.DUP_BRANCHED;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.LINKED;
import static com.hartwig.hmftools.esvee.assembly.RefBaseExtender.extendRefBases;
import static com.hartwig.hmftools.esvee.assembly.AssemblyLinker.tryAssemblyFacing;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.LOCAL_INDEL;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.REMOTE_REGION;
import static com.hartwig.hmftools.esvee.assembly.types.PhaseSet.readsFaceInPhaseSet;
import static com.hartwig.hmftools.esvee.assembly.types.SupportRead.findMatchingFragmentSupport;
import static com.hartwig.hmftools.esvee.assembly.types.SupportRead.hasMatchingFragment;
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
import com.hartwig.hmftools.esvee.assembly.LocalSequenceMatcher;
import com.hartwig.hmftools.esvee.assembly.RemoteRegionAssembler;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.LinkType;
import com.hartwig.hmftools.esvee.assembly.types.PhaseGroup;
import com.hartwig.hmftools.esvee.assembly.types.PhaseSet;
import com.hartwig.hmftools.esvee.assembly.types.RemoteRegion;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;

public class PhaseSetBuilder
{
    private final PhaseGroup mPhaseGroup;
    private final RefGenomeInterface mRefGenome;
    private final RemoteRegionAssembler mRemoteRegionAssembler;
    private final LocalSequenceMatcher mLocalSequenceMatcher;

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
        mLocalSequenceMatcher = new LocalSequenceMatcher(refGenome, LOCAL_ASSEMBLY_MATCH_DISTANCE);

        mPhaseSets = mPhaseGroup.phaseSets();
        mAssemblies = mPhaseGroup.assemblies();
        mSecondarySplitLinks = mPhaseGroup.secondaryLinks();

        mSplitLinks = Lists.newArrayList();
        mFacingLinks = Lists.newArrayList();
    }

    public void buildPhaseSets()
    {
        if(mAssemblies.size() == 1)
        {
            handleSingleAssembly();
            return;
        }

        if(mAssemblies.size() == 2 && handleAssemblyPair())
            return;

        formSplitLinks();

        extendRemoteAssemblies();

        formFacingLinks();

        // extendJunctions();

        formPhaseSets();

        addChainedSupport();

        cleanupAssemblies();
    }

    private boolean handleAssemblyPair()
    {
        // simpler routine without prioritising pairs, facing links or branching - return true if a link is found
        JunctionAssembly assembly1 = mAssemblies.get(0);
        JunctionAssembly assembly2 = mAssemblies.get(1);

        AssemblyLink assemblyLink = checkSplitLink(assembly1, assembly2);

        if(assemblyLink != null)
        {
            buildSplitLink(assembly1, assembly2, false);
            mPhaseSets.add(new PhaseSet(assemblyLink));
            return true;
        }

        return false;
    }

    private void handleSingleAssembly()
    {
        // simpler routine without prioritising pairs, facing links or branching - return true if a link is found
        JunctionAssembly assembly = mAssemblies.get(0);

        if(formsLocalLink(assembly))
            return;

        extendRemoteAssemblies();
    }

    private boolean formsLocalLink(final JunctionAssembly assembly)
    {
        AssemblyLink localRefLink = mLocalSequenceMatcher.tryLocalAssemblyLink(assembly);

        if(localRefLink == null)
            return false;

        assembly.setOutcome(LOCAL_INDEL);

        JunctionAssembly localRefAssembly = localRefLink.otherAssembly(assembly);
        localRefAssembly.setOutcome(LOCAL_INDEL);

        mPhaseGroup.addDerivedAssembly(localRefAssembly);
        mSplitLinks.add(localRefLink);

        // no need to build out the link with matching reads etc since only one assembly has read support
        return true;
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

                int sharedCount = 0;

                for(SupportRead support : assembly1.support())
                {
                    if(hasMatchingFragment(assembly2.support(), support))
                        ++sharedCount;
                }

                // also check candidate reads
                for(SupportRead support : assembly1.candidateSupport())
                {
                    if(hasMatchingFragment(assembly2.support(), support))
                        ++sharedCount;
                }

                if(sharedCount > 0 || isLocalLink)
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

            // test if a link can be made
            if(!alreadyLinked)
            {
                // first check if either assembly matches a local sequence (note both are checked)
                if(!sharedReadPair.IsLocalDelDup)
                {
                    boolean matchesLocal = false;

                    if(formsLocalLink(assembly1))
                    {
                        matchesLocal = true;
                        linkedAssemblies.add(assembly1);
                    }

                    if(formsLocalLink(assembly2))
                    {
                        matchesLocal = true;
                        linkedAssemblies.add(assembly2);
                    }

                    if(matchesLocal)
                        continue;
                }

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
        assembly1.setOutcome(LINKED);
        assembly2.setOutcome(LINKED);

        // look for shared reads between the assemblies, and factor in discordant reads which were only considered candidates until now
        List<SupportRead> matchedCandidates1 = Lists.newArrayList();
        List<SupportRead> matchedCandidates2 = Lists.newArrayList();

        /*
        // list is copied since matched candidates will be removed from repeated matching
        List<SupportRead> candidateSupport2 = Lists.newArrayList(assembly2.candidateSupport());

        // find matching reads, and link reads to each other where possible
        checkMatchingCandidateSupport(assembly2, assembly1.candidateSupport(), candidateSupport2, matchedCandidates1, matchedCandidates2);
        checkMatchingCandidateSupport(assembly1, candidateSupport2, Collections.emptyList(), matchedCandidates2, matchedCandidates1);
        */

        checkMatchingCandidateSupport(assembly2, assembly1.candidateSupport(), assembly2.candidateSupport(), matchedCandidates1, matchedCandidates2);
        checkMatchingCandidateSupport(assembly1, assembly2.candidateSupport(), Collections.emptyList(), matchedCandidates2, matchedCandidates1);

        // build out ref-base assembly support from these non-junction reads - both matched discordant and junction mates
        extendRefBases(assembly1, matchedCandidates1, mRefGenome, allowBranching);
        extendRefBases(assembly2, matchedCandidates2, mRefGenome, allowBranching);
    }

    private static void checkMatchingCandidateSupport(
            final JunctionAssembly otherAssembly,
            final List<SupportRead> candidateSupport, final List<SupportRead> otherCandidateSupport,
            final List<SupportRead> matchedCandidates, final List<SupportRead> otherMatchedCandidates)
    {
        // consider each candidate support read to see if it has a matching read in the other assembly's candidates or junction reads
        int index = 0;
        while(index < candidateSupport.size())
        {
            SupportRead candidateRead = candidateSupport.get(index);

            if(candidateRead.type() == SupportType.JUNCTION_MATE) // added automatically to extend the reference
            {
                candidateSupport.remove(index);
                matchedCandidates.add(candidateRead);
                continue;
            }

            // firsrt check for discordant reads with matching support in the other assembly
            if(hasMatchingFragment(otherAssembly.support(), candidateRead))
            {
                candidateSupport.remove(index);
                matchedCandidates.add(candidateRead);
                continue;
            }

            // then check for candidate & candidate matches
            if(!otherCandidateSupport.isEmpty())
            {

                List<SupportRead> matchedCandidateSupport = findMatchingFragmentSupport(otherCandidateSupport, candidateRead);

                if(!matchedCandidateSupport.isEmpty())
                {
                    candidateSupport.remove(index);
                    matchedCandidates.add(candidateRead);

                    // remove from other's candidates to avoid checking again
                    matchedCandidateSupport.forEach(x -> otherCandidateSupport.remove(x));
                    otherMatchedCandidates.addAll(matchedCandidateSupport);

                    continue;
                }

            }

            ++index;
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
                List<String> localReadIds = assembly.support().stream()
                        .filter(x -> remoteRegion.readIds().contains(x.id()))
                        .map(x -> x.id())
                        .collect(Collectors.toList());

                assembly.candidateSupport().stream()
                        .filter(x -> x.type() == SupportType.CANDIDATE_DISCORDANT)
                        .filter(x -> remoteRegion.readIds().contains(x.id()))
                        .forEach(x -> localReadIds.add(x.id()));

                AssemblyLink assemblyLink = mRemoteRegionAssembler.tryRemoteAssemblyLink(assembly, remoteRegion, localReadIds);

                if(assemblyLink == null)
                    continue;

                JunctionAssembly remoteAssembly = assemblyLink.otherAssembly(assembly);

                remoteAssembly.setOutcome(REMOTE_REGION);

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
                        phaseSet.addAssemblyLinkStart(newLink);
                    else
                        phaseSet.addAssemblyLinkEnd(newLink);

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

    private void addChainedSupport()
    {
        if(mPhaseSets.isEmpty())
            return;

        // look for matched candidate reads spanning proximate breakends and add as support
        for(PhaseSet phaseSet : mPhaseSets)
        {
            if(phaseSet.assemblies().size() <= 2)
                continue;

            for(int i = 0; i < phaseSet.assemblies().size(); ++i)
            {
                JunctionAssembly assembly1 = mAssemblies.get(i);

                List<AssemblyLink> assemblyLinks = phaseSet.findAssemblyLinks(assembly1);

                for(int j = i + 1; j < phaseSet.assemblies().size(); ++j)
                {
                    JunctionAssembly assembly2 = mAssemblies.get(j);

                    // ignore already linked assemblies since their support has been matched, and ignore assemblies in a facing link
                    if(assemblyLinks.stream().anyMatch(x -> x.hasAssembly(assembly2)))
                        continue;

                    addMatchingCandidateSupport(phaseSet, assembly1, assembly2);
                }
            }
        }
    }

    private static void addMatchingCandidateSupport(
            final PhaseSet phaseSet, final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        int assemblyIndex1 = phaseSet.assemblyIndex(assembly1);
        byte assemblyOrientation1 = phaseSet.assemblyOrientation(assembly1);
        int assemblyIndex2 = phaseSet.assemblyIndex(assembly2);
        byte assemblyOrientation2 = phaseSet.assemblyOrientation(assembly2);

        for(int i = 0; i <= 1; ++i)
        {
            final JunctionAssembly assembly = (i == 0) ? assembly1 : assembly2;
            final JunctionAssembly otherAssembly = (i == 0) ? assembly2 : assembly1;

            int index = 0;

            while(index < assembly.candidateSupport().size())
            {
                SupportRead candidateRead = assembly.candidateSupport().get(index);

                if(candidateRead.type() == SupportType.JUNCTION_MATE)
                    continue;

                // only link reads into a supporting fragment across chain links if they face towards each other in the chain
                SupportRead matchedRead = otherAssembly.support().stream()
                        .filter(x -> x.matchesFragment(candidateRead)).findFirst().orElse(null);

                if(readsFaceInPhaseSet(
                        assembly1, assemblyIndex1, assemblyOrientation1, candidateRead,
                        assembly2, assemblyIndex2, assemblyOrientation2, matchedRead))
                {
                    assembly.candidateSupport().remove(index);
                    assembly.addDiscordantSupport(candidateRead);
                    continue;
                }

                // otherwise check for candidate matches
                if(i == 0)
                {
                    SupportRead matchedCandidate = otherAssembly.candidateSupport().stream()
                            .filter(x -> x.matchesFragment(candidateRead)).findFirst().orElse(null);

                    if(readsFaceInPhaseSet(
                            assembly1, assemblyIndex1, assemblyOrientation1, candidateRead,
                            assembly2, assemblyIndex2, assemblyOrientation2, matchedCandidate))
                    {
                        assembly.candidateSupport().remove(index);
                        assembly.addDiscordantSupport(candidateRead);

                        otherAssembly.candidateSupport().remove(matchedCandidate);
                        otherAssembly.addDiscordantSupport(matchedCandidate);

                        continue;
                    }
                }

                ++index;
            }
        }
    }

    private void cleanupAssemblies()
    {
        mAssemblies.forEach(x -> x.clearCandidateSupport()); // no longer required

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
