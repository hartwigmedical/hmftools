package com.hartwig.hmftools.esvee.assembly.phase;

import static java.lang.Math.ceil;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_REF_BASE_MAX_GAP;
import static com.hartwig.hmftools.esvee.AssemblyConstants.LOCAL_ASSEMBLY_MATCH_DISTANCE;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PROXIMATE_DUP_LENGTH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.REMOTE_REGION_REF_MIN_READS;
import static com.hartwig.hmftools.esvee.AssemblyConstants.REMOTE_REGION_REF_MIN_READ_PERCENT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.isLocalAssemblyCandidate;
import static com.hartwig.hmftools.esvee.assembly.RefBaseExtender.checkAddRefBaseRead;
import static com.hartwig.hmftools.esvee.assembly.phase.AssemblyLinker.isAssemblyIndelLink;
import static com.hartwig.hmftools.esvee.assembly.phase.ExtensionType.INDEL;
import static com.hartwig.hmftools.esvee.assembly.phase.ExtensionType.LOCAL_DEL_DUP;
import static com.hartwig.hmftools.esvee.assembly.phase.ExtensionType.REMOTE_REF;
import static com.hartwig.hmftools.esvee.assembly.phase.ExtensionType.SPLIT_LINK;
import static com.hartwig.hmftools.esvee.assembly.phase.ExtensionType.UNMAPPED;
import static com.hartwig.hmftools.esvee.assembly.read.Read.findMatchingFragmentSupport;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyLink.swapAssemblies;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.DUP_BRANCHED;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.LINKED;
import static com.hartwig.hmftools.esvee.assembly.RefBaseExtender.extendRefBases;
import static com.hartwig.hmftools.esvee.assembly.phase.AssemblyLinker.tryAssemblyFacing;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.LOCAL_INDEL;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.REMOTE_LINK;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.REMOTE_REGION;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.SECONDARY;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.UNSET;
import static com.hartwig.hmftools.esvee.assembly.types.SupportRead.hasFragmentOtherRead;
import static com.hartwig.hmftools.esvee.assembly.types.SupportRead.hasMatchingFragmentRead;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.DISCORDANT;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.EXTENSION;
import static com.hartwig.hmftools.esvee.common.CommonUtils.isLineInsertPair;
import static com.hartwig.hmftools.esvee.common.CommonUtils.withLineProximity;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.assembly.AssemblyUtils;
import com.hartwig.hmftools.esvee.assembly.RefBaseExtender;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.LinkType;
import com.hartwig.hmftools.esvee.assembly.types.PhaseGroup;
import com.hartwig.hmftools.esvee.assembly.types.PhaseSet;
import com.hartwig.hmftools.esvee.assembly.types.RemoteRegion;
import com.hartwig.hmftools.esvee.common.CommonUtils;

public class PhaseSetBuilder
{
    private final PhaseGroup mPhaseGroup;
    private final RefGenomeInterface mRefGenome;
    private final RemoteRegionAssembler mRemoteRegionAssembler;
    private final LocalSequenceMatcher mLocalSequenceMatcher;

    // references from phase group
    private final List<JunctionAssembly> mAssemblies;
    private final List<PhaseSet> mPhaseSets; // final proposed phase sets
    private final List<AssemblyLink> mSecondarySplitLinks;
    private final List<ExtensionCandidate> mExtensionCandidates;

    private final Set<JunctionAssembly> mLineRelatedAssemblies;
    private final List<JunctionAssembly> mBranchedAssemblies;

    // working cache only
    private final Set<JunctionAssembly> mLocallyLinkedAssemblies;
    private final List<AssemblyLink> mSplitLinks;
    private final List<AssemblyLink> mFacingLinks;

    // performance tracking
    private double mPerfLogTime;
    private long mStartTimeMs;

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

        mLineRelatedAssemblies = Sets.newHashSet();

        mExtensionCandidates = Lists.newArrayList();
        mSplitLinks = Lists.newArrayList();
        mFacingLinks = Lists.newArrayList();
        mBranchedAssemblies = Lists.newArrayList();
        mLocallyLinkedAssemblies = Sets.newHashSet();

        mPerfLogTime = 0;
        mStartTimeMs = 0;
    }

    public void setPerfLogTime(double perfLogTime) { mPerfLogTime = perfLogTime; }

    public void buildPhaseSets()
    {
        mStartTimeMs = System.currentTimeMillis();

        if(mAssemblies.size() > 100)
        {
            SV_LOGGER.debug("pgId({}) assemblies({}) starting phase set building", mPhaseGroup.id(), mAssemblies.size());
        }

        findLineExtensions();

        findLocalLinks();

        checkLogPerfTime("findLocalLinks");

        findOtherLinksAndExtensions();

        checkLogPerfTime("findOtherLinksAndExtensions");

        addUnlinkedAssemblyRefSupport();

        formFacingLinks();

        checkBranchedAssemblies();

        formPhaseSets();

        addChainedSupport();

        cleanupAssemblies();

        checkLogPerfTime("phaseSets");
    }

    private void findLocalLinks()
    {
        // find local candidate links
        findSplitLinkCandidates(true);

        // prioritise and capture local links
        Collections.sort(mExtensionCandidates, new ExtensionCandidate.LocalLinkComparator());

        for(ExtensionCandidate extensionCandidate : mExtensionCandidates)
        {
            if(!extensionCandidate.isValid())
                continue;

            AssemblyLink assemblyLink = extensionCandidate.Link;

            boolean allowBranching = !(assemblyLink.svType() == DUP && assemblyLink.length() < PROXIMATE_DUP_LENGTH) && mAssemblies.size() > 2;

            applySplitLinkSupport(extensionCandidate.Assembly, extensionCandidate.SecondAssembly, allowBranching);

            extensionCandidate.Assembly.setOutcome(LINKED);
            extensionCandidate.SecondAssembly.setOutcome(LINKED);

            mSplitLinks.add(assemblyLink);
            mLocallyLinkedAssemblies.add(extensionCandidate.Assembly);
            mLocallyLinkedAssemblies.add(extensionCandidate.SecondAssembly);
        }

        // check for other local alignments
        List<JunctionAssembly> unlinkedAssemblies = mAssemblies.stream()
                .filter(x -> !mLocallyLinkedAssemblies.contains(x)).collect(Collectors.toList());

        for(JunctionAssembly assembly : unlinkedAssemblies)
        {
            formsLocalLink(assembly);
        }
    }

    private boolean formsLocalLink(final JunctionAssembly assembly)
    {
        AssemblyLink localRefLink = mLocalSequenceMatcher.tryLocalAssemblyLink(assembly);

        if(localRefLink == null)
            return false;

        assembly.setOutcome(LOCAL_INDEL);

        mLocallyLinkedAssemblies.add(assembly);

        // no need to register this against the phase group, but keep the link since it will be used to create the ref breakend if the aligner doesn't
        // mPhaseGroup.addDerivedAssembly(localRefAssembly);

        mSplitLinks.add(localRefLink);
        JunctionAssembly localRefAssembly = localRefLink.otherAssembly(assembly);
        localRefAssembly.setOutcome(LOCAL_INDEL);
        localRefAssembly.setPhaseGroup(mPhaseGroup);

        return true;
    }

    private void findSplitLinkCandidates(boolean localOnly)
    {
        if(mAssemblies.size() < 2)
            return;

        // test each assembly pair - restricted to local-only links if specified
        // if not, check that a pair hasn't already been tested
        // if no support is found then no need to check the link

        List<ExtensionCandidate> existingCandidates = localOnly ? Collections.emptyList() : Lists.newArrayList(mExtensionCandidates);

        for(int i = 0; i < mAssemblies.size(); ++i)
        {
            JunctionAssembly assembly1 = mAssemblies.get(i);

            // allow linking assemblies to be included, so as to allow secondary links to be found
            // if(mLocallyLinkedAssemblies.contains(assembly1))
            //    continue;

            for(int j = i + 1; j < mAssemblies.size(); ++j)
            {
                JunctionAssembly assembly2 = mAssemblies.get(j);

                // if(mLocallyLinkedAssemblies.contains(assembly2))
                //    continue;

                // avoid a second check of the same pair
                if(existingCandidates.stream().anyMatch(x -> x.Assembly == assembly1 && x.SecondAssembly == assembly2))
                    continue;

                boolean isLocalIndel = isAssemblyIndelLink(assembly1, assembly2);

                boolean isLocalLink = isLocalIndel || isLocalAssemblyCandidate(assembly1, assembly2, false);

                if(localOnly && !isLocalLink)
                    continue;

                boolean hasSharedFragments = hasSharedFragments(assembly1, assembly2);

                AssemblyLink assemblyLink = null;

                if(hasSharedFragments)
                    assemblyLink = checkSplitLink(assembly1, assembly2);

                if(!hasSharedFragments || assemblyLink == null)
                {
                    if(localOnly)
                    {
                        // cache to avoid checking on a second pass
                        mExtensionCandidates.add(new ExtensionCandidate(LOCAL_DEL_DUP, assembly1, assembly2));
                    }

                    continue;
                }

                ExtensionType type = isLocalLink ? (isLocalIndel ? INDEL : LOCAL_DEL_DUP) : SPLIT_LINK;

                ExtensionCandidate extensionCandidate = new ExtensionCandidate(type, assemblyLink);
                mExtensionCandidates.add(extensionCandidate);

                // now count up all possible linking fragments so as to compare with other candidate links and extensions
                Set<String> firstSupportReadIds = assembly1.support().stream().map(x -> x.id()).collect(Collectors.toSet());
                Set<String> firstCandidateReadIds = assembly1.candidateSupport().stream().map(x -> x.id()).collect(Collectors.toSet());

                for(SupportRead support : assembly2.support())
                {
                    if(firstSupportReadIds.contains(support.id()))
                    {
                        firstSupportReadIds.remove(support.id());
                        ++extensionCandidate.AssemblyMatchedSupport;
                        ++extensionCandidate.SecondAssemblyMatchedSupport;
                    }

                    if(firstCandidateReadIds.contains(support.id()))
                    {
                        firstCandidateReadIds.remove(support.id());
                        ++extensionCandidate.AssemblyCandidateReads;
                        ++extensionCandidate.SecondAssemblyMatchedSupport;
                    }
                }

                for(Read read : assembly2.candidateSupport())
                {
                    if(firstSupportReadIds.contains(read.id()))
                    {
                        firstSupportReadIds.remove(read.id());
                        ++extensionCandidate.AssemblyMatchedSupport;
                        ++extensionCandidate.SecondAssemblyCandidateReads;
                    }

                    if(firstCandidateReadIds.contains(read.id()))
                    {
                        firstCandidateReadIds.remove(read.id());
                        ++extensionCandidate.AssemblyCandidateReads;
                        ++extensionCandidate.SecondAssemblyCandidateReads;
                    }
                }
            }
        }
    }

    private void findOtherLinksAndExtensions()
    {
        findUnmappedExtensions();

        findSplitLinkCandidates(false); // since local candidate links have already been found and applied

        findRemoteRefCandidates();

        // prioritise and select from all remaining candidates
        List<ExtensionCandidate> remainingCandidates = mExtensionCandidates.stream()
                .filter(x -> x.isValid())
                .filter(x -> !mLocallyLinkedAssemblies.contains(x.Assembly) && !mLocallyLinkedAssemblies.contains(x.SecondAssembly))
                .collect(Collectors.toList());

        Collections.sort(remainingCandidates, new ExtensionCandidate.StandardComparator());

        Set<JunctionAssembly> primaryLinkedAssemblies = Sets.newHashSet(mLocallyLinkedAssemblies);

        for(ExtensionCandidate extensionCandidate : remainingCandidates)
        {
            if(extensionCandidate.Type == SPLIT_LINK)
            {
                boolean eitherInPrimary = primaryLinkedAssemblies.contains(extensionCandidate.Assembly)
                        || primaryLinkedAssemblies.contains(extensionCandidate.SecondAssembly);

                applySplitLink(extensionCandidate.Link, !eitherInPrimary);

                if(!eitherInPrimary)
                {
                    primaryLinkedAssemblies.add(extensionCandidate.Assembly);
                    primaryLinkedAssemblies.add(extensionCandidate.SecondAssembly);
                }
            }
            else if(extensionCandidate.Type == REMOTE_REF)
            {
                AssemblyLink assemblyLink = extensionCandidate.Link;

                JunctionAssembly initialAssembly = mAssemblies.stream()
                        .filter(x -> x == assemblyLink.first() || x == assemblyLink.second()).findFirst().orElse(null);

                boolean inPrimary = primaryLinkedAssemblies.contains(initialAssembly);

                applyRemoteRefLink(extensionCandidate.Link, initialAssembly, !inPrimary);

                if(!inPrimary)
                    primaryLinkedAssemblies.add(initialAssembly);
            }
            else if(extensionCandidate.Type == UNMAPPED)
            {
                if(!primaryLinkedAssemblies.contains(extensionCandidate.Assembly))
                    applyUnmappedReadExtension(extensionCandidate);
            }
        }
    }

    private void findUnmappedExtensions()
    {
        // any assembly not in a link uses unmapped reads to try to extend the extension sequence
        for(JunctionAssembly assembly : mAssemblies)
        {
            if(assembly.unmappedReads().isEmpty())
                continue;

            if(mLineRelatedAssemblies.contains(assembly)) // ignore if already processed as a line site
                continue;

            UnmappedBaseExtender unmappedBaseExtender = new UnmappedBaseExtender(assembly);
            unmappedBaseExtender.processReads(assembly.unmappedReads());

            if(!unmappedBaseExtender.supportReads().isEmpty())
            {
                ExtensionCandidate extensionCandidate = new ExtensionCandidate(
                        UNMAPPED, assembly, unmappedBaseExtender, unmappedBaseExtender.supportReads().size());

                extensionCandidate.ExtraInfo = format("readSpan(%d)", unmappedBaseExtender.extensionBases().length);
                extensionCandidate.AssemblyCandidateReads = unmappedBaseExtender.supportReads().size();
                mExtensionCandidates.add(extensionCandidate);
            }
        }
    }

    private void findLineExtensions()
    {
        // find line insertion sites - pairs of proximate INDEL-type junctions with a line motif - and jointly search
        // for extension reads across the two assemblies
        mAssemblies.stream().filter(x -> x.hasLineSequence()).forEach(x -> mLineRelatedAssemblies.add(x));

        if(mLineRelatedAssemblies.isEmpty())
            return;

        List<JunctionAssembly> proximateLineAssemblies = Lists.newArrayList();

        Set<JunctionAssembly> extendedAssemblies = Sets.newHashSet();

        for(JunctionAssembly assembly : mAssemblies)
        {
            if(assembly.hasLineSequence())
                continue;

            JunctionAssembly lineAssembly = mLineRelatedAssemblies.stream().filter(x -> isLineInsertPair(assembly, x)).findFirst().orElse(null);

            if(lineAssembly == null)
                continue;

            proximateLineAssemblies.add(assembly);

            List<Read> sharedUnmappedReads = Lists.newArrayList();
            List<RemoteRegion> combinedRemoteRegions = Lists.newArrayList();

            List<JunctionAssembly> lineSiteAssemblies = List.of(assembly, lineAssembly);

            for(JunctionAssembly lineSiteAssembly : lineSiteAssemblies)
            {
                sharedUnmappedReads.addAll(lineSiteAssembly.unmappedReads());

                // collect remote regions if from LINE assemblies or those very close to a LINE assembly
                lineSiteAssembly.remoteRegions().stream()
                        .filter(x -> !x.isSuppOnlyRegion())
                        .forEach(x -> combinedRemoteRegions.add(x));
            }

            RemoteRegion.mergeRegions(combinedRemoteRegions);

            for(RemoteRegion remoteRegion : combinedRemoteRegions)
            {
                List<Read> remoteReads = mRemoteRegionAssembler.extractRemoteReads(remoteRegion);
                sharedUnmappedReads.addAll(remoteReads);
            }

            if(sharedUnmappedReads.isEmpty())
                continue;

            for(JunctionAssembly lineSiteAssembly : lineSiteAssemblies)
            {
                if(extendedAssemblies.contains(lineSiteAssembly))
                    continue;

                extendedAssemblies.add(lineSiteAssembly);

                UnmappedBaseExtender unmappedBaseExtender = new UnmappedBaseExtender(lineSiteAssembly);
                unmappedBaseExtender.processReads(Lists.newArrayList(sharedUnmappedReads)); // list copied so it is given to all assemblies in full

                if(!unmappedBaseExtender.supportReads().isEmpty())
                {
                    lineSiteAssembly.expandExtensionBases(
                            unmappedBaseExtender.extensionBases(), unmappedBaseExtender.baseQualities(), unmappedBaseExtender.supportReads());
                }
            }
        }

        mLineRelatedAssemblies.addAll(proximateLineAssemblies);
    }

    private void findRemoteRefCandidates()
    {
        boolean applyThresholds = mAssemblies.size() > 50;

        for(JunctionAssembly assembly : mAssemblies)
        {
            if(!RemoteRegionAssembler.isExtensionCandidateAssembly(assembly))
                continue;

            // collect remote regions which aren't only supplementaries
            // the check for overlaps with other assemblies in the phase group has been removed since was hiding valid links
            List<RemoteRegion> remoteRegions = assembly.remoteRegions().stream()
                    .filter(x -> !x.isSuppOnlyRegion())
                    .filter(x -> !applyThresholds || x.readIds().size() >= REMOTE_REGION_REF_MIN_READS)
                    .collect(Collectors.toList());

            if(remoteRegions.isEmpty())
                continue;

            // evaluate by remote regions with most linked reads
            Collections.sort(remoteRegions, Comparator.comparingInt(x -> -x.nonSuppReadCount()));

            int minReadCount = applyThresholds ? REMOTE_REGION_REF_MIN_READS : 1;

            if(remoteRegions.size() > 10)
            {
                int maxRemoteReads = remoteRegions.get(0).readCount();
                minReadCount = max(REMOTE_REGION_REF_MIN_READS, (int)ceil(REMOTE_REGION_REF_MIN_READ_PERCENT * maxRemoteReads));
            }

            for(RemoteRegion remoteRegion : remoteRegions)
            {
                if(remoteRegion.readCount() < minReadCount)
                    continue;

                Set<String> localReadIds = assembly.support().stream()
                        .filter(x -> remoteRegion.readIds().contains(x.id()))
                        .map(x -> x.id())
                        .collect(Collectors.toSet());

                int supportCount = localReadIds.size();

                assembly.candidateSupport().stream()
                        .filter(x -> !x.hasJunctionMate())
                        .filter(x -> remoteRegion.hasReadId(x.id()))
                        .forEach(x -> localReadIds.add(x.id()));

                int candidateCount = localReadIds.size() - supportCount;

                if(localReadIds.size() < minReadCount)
                    continue;

                AssemblyLink assemblyLink = mRemoteRegionAssembler.tryRemoteAssemblyLink(assembly, remoteRegion, localReadIds);

                if(assemblyLink == null)
                    continue;

                JunctionAssembly remoteAssembly = assemblyLink.otherAssembly(assembly);

                ExtensionCandidate extensionCandidate = new ExtensionCandidate(REMOTE_REF, assemblyLink);
                extensionCandidate.AssemblyMatchedSupport = supportCount;
                extensionCandidate.AssemblyCandidateReads = candidateCount;
                extensionCandidate.ExtraInfo = format("readSpan(%d)", remoteAssembly.refBaseLength());

                mExtensionCandidates.add(extensionCandidate);
            }
        }
    }

    private void applySplitLink(final AssemblyLink assemblyLink, boolean isPrimaryLink)
    {
        boolean allowBranching = isPrimaryLink && !(assemblyLink.svType() == DUP && assemblyLink.length() < PROXIMATE_DUP_LENGTH);

        applySplitLinkSupport(assemblyLink.first(), assemblyLink.second(), allowBranching);

        if(isPrimaryLink)
        {
            mSplitLinks.add(assemblyLink);

            assemblyLink.first().setOutcome(LINKED);
            assemblyLink.second().setOutcome(LINKED);
        }
        else
        {
            mSecondarySplitLinks.add(assemblyLink);

            // won't override if already set
            assemblyLink.first().setOutcome(SECONDARY);
            assemblyLink.second().setOutcome(SECONDARY);
        }
    }

    private void applyRemoteRefLink(final AssemblyLink assemblyLink, final JunctionAssembly initialAssembly, boolean isPrimaryLink)
    {
        JunctionAssembly remoteAssembly = assemblyLink.otherAssembly(initialAssembly);

        // check for an exact match with an existing assembly, either standard or remote
        JunctionAssembly matchedAssembly = findMatchingAssembly(remoteAssembly, false);

        AssemblyLink remoteLink;

        if(matchedAssembly != null)
        {
            remoteLink = AssemblyLink.swapAssemblies(assemblyLink, remoteAssembly, matchedAssembly);
            remoteAssembly = matchedAssembly;
            isPrimaryLink = false;
        }
        else
        {
            remoteLink = assemblyLink;
            remoteAssembly.setOutcome(REMOTE_REGION);
            mPhaseGroup.addDerivedAssembly(remoteAssembly);
        }

        if(isPrimaryLink)
        {
            // only form one remote link for each assembly
            applySplitLinkSupport(initialAssembly, remoteAssembly, true);
            initialAssembly.setOutcome(REMOTE_LINK);
            mSplitLinks.add(remoteLink);
        }
        else
        {
            applySplitLinkSupport(initialAssembly, remoteAssembly, false);
            mSecondarySplitLinks.add(remoteLink);
        }
    }

    private void applyUnmappedReadExtension(final ExtensionCandidate extensionCandidate)
    {
        UnmappedBaseExtender unmappedBaseExtender = (UnmappedBaseExtender)extensionCandidate.Extender;
        JunctionAssembly assembly = extensionCandidate.Assembly;

        if(!unmappedBaseExtender.supportReads().isEmpty())
        {
            SV_LOGGER.trace("assembly({}) extended {} -> {} with {} unmapped reads",
                    assembly, assembly.extensionLength(), unmappedBaseExtender.extensionBaseLength(),
                    unmappedBaseExtender.supportReads().size());

            assembly.expandExtensionBases(
                    unmappedBaseExtender.extensionBases(), unmappedBaseExtender.baseQualities(), unmappedBaseExtender.supportReads());
        }
    }

    private static boolean hasSharedFragments(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        if(assembly1.support().stream().anyMatch(x -> hasMatchingFragmentRead(assembly2.support(), x)))
            return true;

        if(assembly1.candidateSupport().stream().anyMatch(x -> hasMatchingFragmentRead(assembly2.support(), x)))
            return true;

        if(assembly2.candidateSupport().stream().anyMatch(x -> hasMatchingFragmentRead(assembly1.support(), x)))
            return true;

        if(assembly2.candidateSupport().stream().anyMatch(x -> Read.hasMatchingFragmentRead(assembly1.candidateSupport(), x)))
            return true;

        return false;
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

    private void addUnlinkedAssemblyRefSupport()
    {
        // any assembly which did not form a link or only an unmapped extension will now extend its ref bases from junction & extension mates
        for(JunctionAssembly assembly : mAssemblies)
        {
            if(assembly.outcome() != UNSET)
                continue;
            
            if(mSplitLinks.stream().anyMatch(x -> x.hasAssembly(assembly)))
                continue;

            // add junction mate reads to the ref side bases
            List<Read> refExtensionReads = Lists.newArrayList();

            List<SupportRead> extensionSupport = assembly.support().stream().filter(x -> x.type() == EXTENSION).collect(Collectors.toList());

            for(Read read : assembly.candidateSupport())
            {
                if(read.hasJunctionMate())
                {
                    refExtensionReads.add(read);
                }
                else
                {
                    if(extensionSupport.stream().anyMatch(x -> x.matchesFragment(read, false)))
                        refExtensionReads.add(read);
                }
            }

            extendRefBases(assembly, refExtensionReads, mRefGenome, false, false);
        }
    }

    private boolean applySplitLinkSupport(final JunctionAssembly assembly1, final JunctionAssembly assembly2, boolean allowBranching)
    {
        // look for shared reads between the assemblies, and factor in discordant reads which were only considered candidates until now
        List<Read> matchedCandidates1 = Lists.newArrayList();
        List<Read> matchedCandidates2 = Lists.newArrayList();

        addLocalMateSupport(assembly1, assembly2);

        checkMatchingCandidateSupport(assembly2, assembly1.candidateSupport(), assembly2.candidateSupport(), matchedCandidates1, matchedCandidates2);
        checkMatchingCandidateSupport(assembly1, assembly2.candidateSupport(), Collections.emptyList(), matchedCandidates2, matchedCandidates1);

        addMatchingExtensionCandidates(assembly1, matchedCandidates1);
        addMatchingExtensionCandidates(assembly2, matchedCandidates2);

        // remove any ref discordant candidates if their only criteria for inclusion is being long
        List<Read> refCandidates1 = Lists.newArrayList();
        boolean hasNonLocalTumorFragment = false;
        boolean hasNonLocalRefFragment = false;

        for(Read read : matchedCandidates1)
        {
            if(read.isReference())
            {
                refCandidates1.add(read);

                hasNonLocalRefFragment |= CommonUtils.isDiscordantFragment(
                        read.bamRecord(), -1, read.supplementaryData());
            }
            else
            {
                hasNonLocalTumorFragment |= CommonUtils.isDiscordantFragment(
                        read.bamRecord(), -1, read.supplementaryData());
            }
        }

        if(hasNonLocalTumorFragment && !hasNonLocalRefFragment)
        {
            List<Read> refCandidates2 = matchedCandidates2.stream().filter(x -> x.isReference()).collect(Collectors.toList());
            refCandidates1.forEach(x -> matchedCandidates1.remove(x));
            refCandidates2.forEach(x -> matchedCandidates2.remove(x));
        }

        if(matchedCandidates1.isEmpty() && matchedCandidates2.isEmpty())
            return false;

        // build out ref-base assembly support from these non-junction reads - both matched discordant and junction mates
        extendRefBases(assembly1, matchedCandidates1, mRefGenome, allowBranching, true);
        extendRefBases(assembly2, matchedCandidates2, mRefGenome, allowBranching, true);

        for(JunctionAssembly assembly : mAssemblies)
        {
            if(assembly.outcome() == DUP_BRANCHED && !mBranchedAssemblies.contains(assembly))
                mBranchedAssemblies.add(assembly);
        }

        return true;
    }

    private static void checkMatchingCandidateSupport(
            final JunctionAssembly otherAssembly,
            final List<Read> candidateSupport, final List<Read> otherCandidateSupport,
            final List<Read> matchedCandidates, final List<Read> otherMatchedCandidates)
    {
        // consider each candidate support read to see if it has a matching read in the other assembly's candidates or junction reads
        int index = 0;
        while(index < candidateSupport.size())
        {
            Read candidateRead = candidateSupport.get(index);

            if(candidateRead.hasJunctionMate()) // added automatically to extend the reference
            {
                candidateSupport.remove(index);
                matchedCandidates.add(candidateRead);
                continue;
            }

            // first check for discordant reads with matching support in the other assembly
            if(hasFragmentOtherRead(otherAssembly.support(), candidateRead))
            {
                candidateSupport.remove(index);
                matchedCandidates.add(candidateRead);
                continue;
            }

            // then check for candidate & candidate matches
            if(!otherCandidateSupport.isEmpty())
            {
                List<Read> matchedCandidateSupport = findMatchingFragmentSupport(otherCandidateSupport, candidateRead);

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

    private static void addMatchingExtensionCandidates(final JunctionAssembly assembly, final List<Read> matchedCandidates)
    {
        Set<String> extensionSupport = assembly.support().stream()
                .filter(x -> x.type() == EXTENSION).map(x -> x.id()).collect(Collectors.toSet());

        if(extensionSupport.isEmpty())
            return;

        int index = 0;
        while(index < assembly.candidateSupport().size())
        {
            Read candidateRead = assembly.candidateSupport().get(index);

            if(extensionSupport.remove(candidateRead.id()))
            {
                matchedCandidates.add(candidateRead);
                assembly.candidateSupport().remove(index);
                continue;
            }

            ++index;
        }
    }

    private void addLocalMateSupport(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        if(!isLocalAssemblyCandidate(assembly1, assembly2, false))
            return;

        // look for concordant mate reads which are on the otherr side of the junction and so were initially excluded
        for(int i = 0; i <= 1; ++i)
        {
            JunctionAssembly assembly = (i == 0) ? assembly1 : assembly2;
            JunctionAssembly otherAssembly = (i == 0) ? assembly2 : assembly1;

            for(Read mateRead : assembly.concordantCandidates())
            {
                // check mate orientation and position vs the other assembly's junction
                if(otherAssembly.isForwardJunction())
                {
                    if(mateRead.orientation().isForward() && mateRead.alignmentEnd() <= otherAssembly.junction().Position)
                    {
                        otherAssembly.addCandidateSupport(mateRead);
                    }
                }
                else
                {
                    if(mateRead.orientation().isReverse() && mateRead.alignmentStart() >= otherAssembly.junction().Position)
                    {
                        otherAssembly.addCandidateSupport(mateRead);
                    }
                }
            }
        }
    }

    private void formFacingLinks()
    {
        if(mAssemblies.size() == 1 || (mAssemblies.size() == 2 && mSplitLinks.size() == 1))
            return;

        // for each assembly in a split link, look for a facing link (whether linked or not)
        Set<JunctionAssembly> facingAssemblies = Sets.newHashSet();

        for(int i = 0; i < mAssemblies.size() - 1; ++i)
        {
            JunctionAssembly assembly1 = mAssemblies.get(i);

            for(int j = i + 1; j < mAssemblies.size(); ++j)
            {
                JunctionAssembly assembly2 = mAssemblies.get(j);

                if(facingAssemblies.contains(assembly1) || facingAssemblies.contains(assembly2))
                    continue;

                AssemblyLink facingLink = tryAssemblyFacing(assembly1, assembly2);

                if(facingLink == null)
                    continue;

                // compelling evidence is a read from the new assembly which overlaps with the linked junction's reads
                // if(assembliesShareReads(assembly2, splitAssembly))
                mFacingLinks.add(facingLink);
                facingAssemblies.add(assembly1);
                facingAssemblies.add(assembly2);
            }
        }
    }

    private void checkBranchedAssemblies()
    {
        for(JunctionAssembly assembly : mBranchedAssemblies)
        {
            boolean hasSplitLink = mSplitLinks.stream().anyMatch(x -> x.hasAssembly(assembly));
            boolean hasFacingLink = mFacingLinks.stream().anyMatch(x -> x.hasAssembly(assembly));

            if(!hasSplitLink && !hasFacingLink) // will be removed later on
                continue;

            // check if the original assembly could be swapped out for this one if it's in a facing link but not a split link
            if(hasFacingLink && !hasSplitLink)
            {
                JunctionAssembly originalAssembly = findMatchingAssembly(assembly, true);

                if(originalAssembly == null)
                    continue;

                AssemblyLink existingLink = mSplitLinks.stream().filter(x -> x.hasAssembly(originalAssembly)).findFirst().orElse(null);
                boolean origHasFacingLink = mFacingLinks.stream().anyMatch(x -> x.hasAssembly(originalAssembly));

                if(!origHasFacingLink && existingLink != null)
                {
                    // swap the assemblies
                    AssemblyLink replacementLink = swapAssemblies(existingLink, originalAssembly, assembly);
                    mSplitLinks.remove(existingLink);
                    mSplitLinks.add(replacementLink);
                    mPhaseGroup.derivedAssemblies().remove(assembly);
                    originalAssembly.setOutcome(DUP_BRANCHED, true);
                    assembly.setOutcome(LINKED);
                }
            }
        }
    }

    private JunctionAssembly findMatchingAssembly(final JunctionAssembly assembly, boolean requireExtensionMatch)
    {
        return AssemblyUtils.findMatchingAssembly(mAssemblies, assembly, requireExtensionMatch);
    }

    private void formPhaseSets()
    {
        // use split and facing links to assign assemblies to phase sets
        while(!mSplitLinks.isEmpty() || !mFacingLinks.isEmpty())
        {
            AssemblyLink assemblyLink = !mSplitLinks.isEmpty() ? mSplitLinks.remove(0) : mFacingLinks.remove(0);

            PhaseSet phaseSet = new PhaseSet(assemblyLink);
            mPhaseSets.add(phaseSet);

            // look for facing and then splits links for this phase set
            for(int se = SE_START; se <= SE_END; ++se)
            {
                // check start and then end links of this phase set
                JunctionAssembly linkingAssembly = (se == SE_START) ? assemblyLink.first() : assemblyLink.second();
                boolean findSplit = assemblyLink.type() == LinkType.FACING;

                while(true)
                {
                    AssemblyLink nextLink = findLinkedAssembly(linkingAssembly, findSplit);

                    if(nextLink == null)
                        break;

                    if(se == SE_START)
                        phaseSet.addAssemblyLinkStart(nextLink);
                    else
                        phaseSet.addAssemblyLinkEnd(nextLink);

                    findSplit = !findSplit;
                    linkingAssembly = nextLink.otherAssembly(linkingAssembly);
                }
            }
        }

        Set<JunctionAssembly> processedSecondaries = Sets.newHashSet();

        for(AssemblyLink link : mSecondarySplitLinks)
        {
            boolean firstIsLineSite = link.first().outcome() == SECONDARY && mLineRelatedAssemblies.contains(link.first());
            boolean secondIsLineSite = link.second().outcome() == SECONDARY && mLineRelatedAssemblies.contains(link.second());

            if(firstIsLineSite || secondIsLineSite)
            {
                boolean processed = (firstIsLineSite && processedSecondaries.contains(link.first()))
                        || (secondIsLineSite && processedSecondaries.contains(link.second()));

                if(!processed)
                {
                    PhaseSet phaseSet = new PhaseSet(link);
                    mPhaseSets.add(phaseSet);

                    if(firstIsLineSite)
                        processedSecondaries.add(link.first());

                    if(secondIsLineSite)
                        processedSecondaries.add(link.second());
                }
            }
            else
            {
                // add to the the related phase
                for(PhaseSet phaseSet : mPhaseSets)
                {
                    if(phaseSet.hasAssembly(link.first()) || phaseSet.hasAssembly(link.second()))
                        phaseSet.addSecondaryLink(link);
                }
            }
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

            for(int i = 0; i < phaseSet.assemblies().size() - 1; ++i)
            {
                JunctionAssembly assembly1 = phaseSet.assemblies().get(i);

                List<AssemblyLink> assemblyLinks = phaseSet.findAssemblyLinks(assembly1);

                for(int j = i + 1; j < phaseSet.assemblies().size(); ++j)
                {
                    JunctionAssembly assembly2 = phaseSet.assemblies().get(j);

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
        // assemblies must face each other in the chain
        if(!phaseSet.assembliesFaceInPhaseSet(assembly1, assembly2))
            return;

        for(int i = 0; i <= 1; ++i)
        {
            final JunctionAssembly assembly = (i == 0) ? assembly1 : assembly2;
            final JunctionAssembly otherAssembly = (i == 0) ? assembly2 : assembly1;

            int index = 0;

            while(index < assembly.candidateSupport().size())
            {
                Read candidateRead = assembly.candidateSupport().get(index);

                if(candidateRead.hasJunctionMate())
                {
                    ++index;
                    continue;
                }

                // only link reads into a supporting fragment across chain links if they face towards each other in the chain
                SupportRead matchedRead = otherAssembly.support().stream()
                        .filter(x -> x.matchesFragment(candidateRead, false)).findFirst().orElse(null);

                if(matchedRead != null)
                {
                    assembly.candidateSupport().remove(index);
                    checkAddRefBaseRead(assembly, candidateRead, DISCORDANT);
                    continue;
                }

                // otherwise check for candidate matches
                if(i == 0)
                {
                    Read matchedCandidate = otherAssembly.candidateSupport().stream()
                            .filter(x -> x.matchesFragment(candidateRead, false)).findFirst().orElse(null);

                    if(matchedCandidate != null)
                    {
                        assembly.candidateSupport().remove(index);
                        checkAddRefBaseRead(assembly, candidateRead, DISCORDANT);

                        otherAssembly.candidateSupport().remove(matchedCandidate);
                        checkAddRefBaseRead(otherAssembly, matchedCandidate, DISCORDANT);

                        continue;
                    }
                }

                ++index;
            }
        }
    }

    private void cleanupAssemblies()
    {
        List<JunctionAssembly> branchedAssembliesToRemove = null;

        for(JunctionAssembly assembly : mAssemblies)
        {
            assembly.clearCandidateSupport(); // no further use for candidate reads

            assembly.clearSupportCachedReads(); // remove references to actual SAMRecords, keeping only summary info

            boolean inPhaseSet = mPhaseSets.stream().anyMatch(x -> x.hasAssembly(assembly));

            if(inPhaseSet && assembly.outcome() == UNSET)
                assembly.setOutcome(LINKED);

            RefBaseExtender.trimAssemblyRefBases(assembly, ASSEMBLY_REF_BASE_MAX_GAP);

            if(assembly.outcome() == DUP_BRANCHED)
            {
                // remove any branched assemblies which did not form a facing link
                boolean inFacingLink = false;

                if(inPhaseSet)
                {
                    for(PhaseSet phaseSet : mPhaseSets)
                    {
                        if(phaseSet.assemblyLinks().stream().filter(x -> x.type() == LinkType.FACING).anyMatch(x -> x.hasAssembly(assembly)))
                        {
                            // set outcome to original assembly
                            JunctionAssembly originalAssembly = findMatchingAssembly(assembly, false);

                            if(originalAssembly != null)
                                assembly.setOutcome(originalAssembly.outcome());

                            inFacingLink = true;
                            break;
                        }
                    }
                }

                if(!inFacingLink)
                {
                    if(branchedAssembliesToRemove == null)
                        branchedAssembliesToRemove = Lists.newArrayList(assembly);
                    else
                        branchedAssembliesToRemove.add(assembly);
                }
            }
        }

        // finally remove any branched assemblies which did not form a facing link
        if(branchedAssembliesToRemove != null)
        {
            for(JunctionAssembly branchedAssembly : branchedAssembliesToRemove)
            {
                mPhaseGroup.assemblies().remove(branchedAssembly);
                mPhaseGroup.derivedAssemblies().remove(branchedAssembly);
            }
        }
    }

    private void checkLogPerfTime(final String stage)
    {
        if(mPerfLogTime == 0)
            return;

        long timeTakenMs = System.currentTimeMillis() - mStartTimeMs;
        double seconds = timeTakenMs / 1000.0;

        if(seconds >= mPerfLogTime)
        {
            StringJoiner sj = new StringJoiner(";");
            for(int i = 0; i < min(mAssemblies.size(), 4); ++i)
            {
                sj.add(mAssemblies.get(i).junction().coords());
            }

            SV_LOGGER.debug(format("pgId(%d) assemblies(%d: %s) stage(%s) time(%.3fs) details(links=%d candidates=%d line=%d remoteRefReads=%d)",
                    mPhaseGroup.id(), mAssemblies.size(), sj, stage, seconds, mSplitLinks.size(), mExtensionCandidates.size(),
                    mLineRelatedAssemblies.size(), mRemoteRegionAssembler != null ? mRemoteRegionAssembler.totalRemoteReadsSearch() : 0));
        }

        mStartTimeMs = System.currentTimeMillis();
    }
}
