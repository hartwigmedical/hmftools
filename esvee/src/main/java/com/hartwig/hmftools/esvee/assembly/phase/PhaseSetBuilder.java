package com.hartwig.hmftools.esvee.assembly.phase;

import static java.lang.Math.ceil;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_REF_BASE_MAX_GAP;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.LOCAL_ASSEMBLY_MATCH_DISTANCE;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PROXIMATE_DUP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.REMOTE_REGION_REF_MIN_READS;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.REMOTE_REGION_REF_MIN_READ_PERCENT;
import static com.hartwig.hmftools.esvee.assembly.RefBaseExtender.checkAddRefBaseRead;
import static com.hartwig.hmftools.esvee.assembly.phase.AssemblyLinker.isAssemblyIndelLink;
import static com.hartwig.hmftools.esvee.assembly.phase.AssemblyLinker.isFacingAssemblyCandidate;
import static com.hartwig.hmftools.esvee.assembly.phase.ExtensionType.INDEL;
import static com.hartwig.hmftools.esvee.assembly.phase.ExtensionType.LOCAL_DEL_DUP;
import static com.hartwig.hmftools.esvee.assembly.phase.ExtensionType.REMOTE_REF;
import static com.hartwig.hmftools.esvee.assembly.phase.ExtensionType.SPLIT_LINK;
import static com.hartwig.hmftools.esvee.assembly.phase.ExtensionType.UNMAPPED;
import static com.hartwig.hmftools.esvee.assembly.phase.RemoteRegionAssembler.collectCandidateRemoteRegions;
import static com.hartwig.hmftools.esvee.assembly.read.Read.findMatchingFragmentSupport;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.isDiscordantFragment;
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
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.DISCORDANT;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.EXTENSION;
import static com.hartwig.hmftools.esvee.common.CommonUtils.isLineInsertPair;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.assembly.AssemblyConfig;
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

    private static final int HIGH_ASSEMBLY_COUNT = 100;
    private static final int HIGH_ASSEMBLY_READ_COUNT = 100;

    public void setPerfLogTime(double perfLogTime) { mPerfLogTime = perfLogTime; }

    public void buildPhaseSets()
    {
        mStartTimeMs = System.currentTimeMillis();

        if(hasHighAssemblyCount())
        {
            SV_LOGGER.debug("pgId({}) assemblies({}) starting phase set building", mPhaseGroup.id(), mAssemblies.size());
        }

        findLineExtensions();

        findLocalLinks();

        findOtherLinksAndExtensions();

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
            boolean allowDiscordantReads = !CommonUtils.isShortLocalDelDupIns(assemblyLink.svType(), assemblyLink.length());

            extensionCandidate.markSelected();
            applySplitLinkSupport(extensionCandidate.Assembly, extensionCandidate.SecondAssembly, allowBranching, allowDiscordantReads);

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

        checkLogPerfTime("findLocalLinks");
    }

    private boolean formsLocalLink(final JunctionAssembly assembly)
    {
        if(assembly.discordantOnly())
            return false;

        // look for a local sequence match for the extension bases, thereby forming a short DEL or DUP
        AssemblyLink localRefLink = mLocalSequenceMatcher.tryLocalAssemblyLink(assembly);

        if(localRefLink == null)
            return false;

        assembly.setOutcome(LOCAL_INDEL);

        mLocallyLinkedAssemblies.add(assembly);

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

            Set<String> firstReadIds = Sets.newHashSet();
            boolean firstHighReadCount = assemblyHasHighReadCount(assembly1);
            populateReadIds(assembly1, firstReadIds, localOnly);

            // allow linking assemblies to be included, to allow secondary links to be found
            for(int j = i + 1; j < mAssemblies.size(); ++j)
            {
                JunctionAssembly assembly2 = mAssemblies.get(j);

                // avoid a second check of the same pair
                if(existingCandidates.stream().anyMatch(x -> x.matchesAssemblies(assembly1, assembly2)))
                    continue;

                boolean isLocalIndel = false;
                boolean isLocalLink = false;

                if(localOnly)
                {
                    isLocalIndel = isAssemblyIndelLink(assembly1, assembly2);

                    isLocalLink = isLocalIndel || isLocalAssemblyCandidate(assembly1, assembly2);

                    if(!isLocalLink)
                        continue;

                    // discordant pairs cannot by definition be short local indel-type links
                    if(assembly1.discordantOnly() || assembly2.discordantOnly())
                        continue;

                    // only link indel junctions to each other
                    if(assembly1.indel() != assembly2.indel())
                        continue;
                }

                Set<String> secondReadIds = Sets.newHashSet();
                boolean secondHighReadCount = assemblyHasHighReadCount(assembly2);
                populateReadIds(assembly2, secondReadIds, localOnly);

                Set<String> firstReadIdsRelated;
                Set<String> secondReadIdsRelated;

                if(!localOnly && (firstHighReadCount || secondHighReadCount))
                {
                    firstReadIdsRelated = Sets.newHashSet();
                    secondReadIdsRelated = Sets.newHashSet();
                    populateNonLocalReadIds(assembly1, assembly2, firstReadIdsRelated, secondReadIdsRelated);
                }
                else
                {
                    firstReadIdsRelated = firstReadIds;
                    secondReadIdsRelated = secondReadIds;
                }

                // proximate breakends may not share reads esp if indels vs soft-clips are the source of differences
                boolean hasSharedFragments = localOnly || hasSharedFragment(firstReadIdsRelated, secondReadIdsRelated);

                AssemblyLink assemblyLink = null;

                if(hasSharedFragments)
                    assemblyLink = checkSplitLink(assembly1, assembly2, isLocalLink);

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

                // now count up all possible linking fragments to compare with other candidate links and extensions
                countSharedFragments(extensionCandidate, firstReadIdsRelated, secondReadIdsRelated);
            }
        }
    }

    private static boolean assemblyHasHighReadCount(final JunctionAssembly assembly)
    {
        return assembly.supportCount() >= HIGH_ASSEMBLY_READ_COUNT || assembly.candidateSupport().size() >= HIGH_ASSEMBLY_READ_COUNT;
    }

    private static void populateReadIds(final JunctionAssembly assembly, final Set<String> readIds, boolean localOnly)
    {
        if(assemblyHasHighReadCount(assembly))
        {
            if(localOnly)
            {
                assembly.support().stream().filter(x -> !x.isDiscordant() && x.isMateMapped()).forEach(x -> readIds.add(x.id()));
                assembly.candidateSupport().stream().filter(x -> !isDiscordantFragment(x) && x.isMateMapped()).forEach(x -> readIds.add(x.id()));
            }
            else
            {
                return; // will be populated based on another non-local assembly
            }
        }
        else
        {
            // take all reads
            assembly.support().forEach(x -> readIds.add(x.id()));
            assembly.candidateSupport().forEach(x -> readIds.add(x.id()));
        }
    }

    private static void populateNonLocalReadIds(
            final JunctionAssembly firstAssembly, final JunctionAssembly secondAssembly,
            final Set<String> firstReadIds, final Set<String> secondReadIds)
    {
        // collect read IDs only from remote regions which overlap each other
        List<RemoteRegion> firstMatchedRegions = firstAssembly.remoteRegions().stream()
                .filter(x -> x.overlapsAssembly(secondAssembly)).collect(Collectors.toList());

        List<RemoteRegion> secondMatchedRegions = secondAssembly.remoteRegions().stream()
                .filter(x -> x.overlapsAssembly(firstAssembly)).collect(Collectors.toList());

        if(firstMatchedRegions.isEmpty() || secondMatchedRegions.isEmpty())
            return;

        firstMatchedRegions.forEach(x -> firstReadIds.addAll(x.readIds()));
        secondMatchedRegions.forEach(x -> secondReadIds.addAll(x.readIds()));
    }

    private static boolean hasSharedFragment(final Set<String> firstReadIds, final Set<String> secondReadIds)
    {
        return firstReadIds.stream().anyMatch(x -> secondReadIds.contains(x));
    }

    private void countSharedFragments(
            final ExtensionCandidate extensionCandidate, final Set<String> firstReadIds, final Set<String> secondReadIds)
    {
        for(String readId : secondReadIds)
        {
            if(firstReadIds.contains(readId))
            {
                ++extensionCandidate.SupportCount;
            }
        }
    }

    private void findOtherLinksAndExtensions()
    {
        findUnmappedExtensions();

        checkLogPerfTime("findUnmappedExtensions");

        findSplitLinkCandidates(false); // since local candidate links have already been found and applied

        findRemoteRefCandidates();

        // prioritise and select from all remaining candidates
        List<ExtensionCandidate> remainingCandidates = mExtensionCandidates.stream()
                .filter(x -> !x.selected()) // skip those already registered
                .filter(x -> x.isValid())
                .collect(Collectors.toList());

        if(remainingCandidates.isEmpty())
            return;

        Collections.sort(remainingCandidates, new ExtensionCandidate.StandardComparator());

        Set<JunctionAssembly> primaryLinkedAssemblies = Sets.newHashSet(mLocallyLinkedAssemblies);

        for(ExtensionCandidate extensionCandidate : remainingCandidates)
        {
            if(extensionCandidate.Type == SPLIT_LINK)
            {
                boolean eitherInPrimary = primaryLinkedAssemblies.contains(extensionCandidate.Assembly)
                        || primaryLinkedAssemblies.contains(extensionCandidate.SecondAssembly);

                if(eitherInPrimary
                && (extensionCandidate.Assembly.discordantOnly() || extensionCandidate.SecondAssembly.discordantOnly()))
                {
                    // ignore discordant-only junctions with links to an assembly already link
                    continue;
                }

                extensionCandidate.markSelected();
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

                extensionCandidate.markSelected();
                applyRemoteRefLink(extensionCandidate.Link, initialAssembly, !inPrimary);

                if(!inPrimary)
                    primaryLinkedAssemblies.add(initialAssembly);
            }
            else if(extensionCandidate.Type == UNMAPPED)
            {
                // extend the assembly if not already linked in any way - this facilitates further links and/or a better SGL for alignment
                if(extensionCandidate.Assembly.outcome() == UNSET)
                {
                    extensionCandidate.markSelected();
                    applyUnmappedReadExtension(extensionCandidate);
                }
            }
        }

        checkLogPerfTime("findOtherLinks");
    }

    private void findUnmappedExtensions()
    {
        // any assembly not in a link uses unmapped reads to try to extend the extension sequence
        boolean hasHighAssemblyCount = hasHighAssemblyCount();

        for(JunctionAssembly assembly : mAssemblies)
        {
            if(mLineRelatedAssemblies.contains(assembly)) // ignore if already processed as a line site
                continue;

            List<Read> unmappedReads = Lists.newArrayList(assembly.unmappedReads());

            if(!AssemblyConfig.RunRemoteRefLinking)
            {
                List<RemoteRegion> combinedRemoteRegions = collectCandidateRemoteRegions(assembly, mAssemblies, hasHighAssemblyCount);
                mRemoteRegionAssembler.extractRemoteRegionReads(mPhaseGroup.id(), combinedRemoteRegions, unmappedReads, hasHighAssemblyCount);
            }

            if(unmappedReads.isEmpty())
                continue;

            UnmappedBaseExtender unmappedBaseExtender = new UnmappedBaseExtender(assembly);
            unmappedBaseExtender.processReads(unmappedReads);

            if(!unmappedBaseExtender.supportReads().isEmpty())
            {
                ExtensionCandidate extensionCandidate = new ExtensionCandidate(
                        UNMAPPED, assembly, unmappedBaseExtender, unmappedBaseExtender.supportReads().size());

                extensionCandidate.ExtraInfo = format("readSpan(%d)", unmappedBaseExtender.extensionBases().length);
                extensionCandidate.SupportCount = unmappedBaseExtender.supportReads().size();
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

            mRemoteRegionAssembler.extractRemoteRegionReads(
                    mPhaseGroup.id(), combinedRemoteRegions, sharedUnmappedReads, hasHighAssemblyCount());

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

        checkLogPerfTime("findLineExtensions");
    }

    private void applySplitLink(final AssemblyLink assemblyLink, boolean isPrimaryLink)
    {
        boolean allowBranching = isPrimaryLink && !(assemblyLink.svType() == DUP && assemblyLink.length() < PROXIMATE_DUP_LENGTH);
        boolean allowDiscordantReads = !CommonUtils.isShortLocalDelDupIns(assemblyLink.svType(), assemblyLink.length());

        applySplitLinkSupport(assemblyLink.first(), assemblyLink.second(), allowBranching, allowDiscordantReads);

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

    private AssemblyLink checkSplitLink(final JunctionAssembly assembly1, final JunctionAssembly assembly2, boolean isLocalIndel)
    {
        if(assembly1.junction() == assembly2.junction()) // ignore duplicates
            return null;

        if(isLocalIndel)
        {
            // handle local INDELs here since the following logic currently applies to them
            AssemblyLink assemblyLink = AssemblyLinker.tryAssemblyIndel(assembly1, assembly2);

            if(assemblyLink != null)
                return assemblyLink;
        }

        return AssemblyLinker.tryAssemblyOverlap(assembly1, assembly2, true, isLocalIndel);
    }

    private void addUnlinkedAssemblyRefSupport()
    {
        // any assembly which did not form a link or only an unmapped extension will now extend its ref bases from junction & extension mates
        List<JunctionAssembly> assemblies = Lists.newArrayList(mAssemblies);

        boolean allowRefSideSoftClipBranching = !AssemblyConfig.RunRemoteRefLinking;

        for(JunctionAssembly assembly : assemblies)
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

            extendRefBases(assembly, refExtensionReads, mRefGenome, allowRefSideSoftClipBranching, allowRefSideSoftClipBranching);
        }
    }

    private boolean applySplitLinkSupport(
            final JunctionAssembly assembly1, final JunctionAssembly assembly2, boolean allowBranching, boolean allowDiscordantReads)
    {
        // look for shared reads between the assemblies, and factor in discordant reads which were only considered candidates until now
        List<Read> matchedCandidates1 = Lists.newArrayList();
        List<Read> matchedCandidates2 = Lists.newArrayList();

        addLocalMateSupport(assembly1, assembly2);

        checkMatchingCandidateSupport(
                assembly2, allowDiscordantReads, assembly1.candidateSupport(), assembly2.candidateSupport(), matchedCandidates1, matchedCandidates2);

        checkMatchingCandidateSupport(
                assembly1, allowDiscordantReads, assembly2.candidateSupport(), Collections.emptyList(), matchedCandidates2, matchedCandidates1);

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

        // register any newly branched assemblies
        for(JunctionAssembly assembly : mAssemblies)
        {
            if(assembly.outcome() == DUP_BRANCHED && !mBranchedAssemblies.contains(assembly))
                mBranchedAssemblies.add(assembly);
        }

        return true;
    }

    private static void checkMatchingCandidateSupport(
            final JunctionAssembly otherAssembly, boolean allowDiscordantReads,
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

            if(!allowDiscordantReads)
            {
                ++index;
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

    protected static boolean isLocalAssemblyCandidate(final JunctionAssembly first, final JunctionAssembly second)
    {
        // just checks orientation, no read concordance
        return AssemblyUtils.isLocalAssemblyCandidate(first, second, false, false);
    }

    private void addLocalMateSupport(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        if(!isLocalAssemblyCandidate(assembly1, assembly2))
            return;

        // look for concordant mate reads which are on the other side of the junction and so were initially excluded
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

            if(!isFacingAssemblyCandidate(assembly1, facingAssemblies, mSplitLinks))
                continue;

            for(int j = i + 1; j < mAssemblies.size(); ++j)
            {
                JunctionAssembly assembly2 = mAssemblies.get(j);

                if(!isFacingAssemblyCandidate(assembly2, facingAssemblies, mSplitLinks))
                    continue;

                AssemblyLink facingLink = tryAssemblyFacing(assembly1, assembly2, mSplitLinks);

                if(facingLink == null)
                    continue;

                // compelling evidence is a read from the new assembly which overlaps with the linked junction's reads
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

            // keep local ref links separate from other assemblies, and these won't be aligned if short
            if(assemblyLink.type() == LinkType.SPLIT && assemblyLink.first().outcome() == LOCAL_INDEL)
                continue;

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

        List<PhaseSet> secondaryLinePhaseSets = Lists.newArrayList();

        for(AssemblyLink link : mSecondarySplitLinks)
        {
            boolean firstIsLineSite = link.first().outcome() == SECONDARY && mLineRelatedAssemblies.contains(link.first());
            boolean secondIsLineSite = link.second().outcome() == SECONDARY && mLineRelatedAssemblies.contains(link.second());

            // put LINE links not already associated with a phase set into their own
            if(firstIsLineSite || secondIsLineSite)
            {
                boolean processed = (firstIsLineSite && processedSecondaries.contains(link.first()))
                        || (secondIsLineSite && processedSecondaries.contains(link.second()));

                if(!processed)
                {
                    PhaseSet phaseSet = new PhaseSet(link);
                    secondaryLinePhaseSets.add(phaseSet);

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

        mPhaseSets.addAll(secondaryLinePhaseSets);
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
            if(!AssemblyConfig.WriteCandidateReads)
                assembly.clearCandidateSupport(); // no further use for candidate reads

            assembly.clearSupportCachedReads(); // remove references to actual SAMRecords, keeping only summary info

            boolean inPhaseSet = false;
            boolean inFacingLink = false;

            for(PhaseSet phaseSet : mPhaseSets)
            {
                if(phaseSet.hasAssembly(assembly))
                {
                    inPhaseSet = true;
                    inFacingLink |= phaseSet.assemblyLinks().stream().anyMatch(x -> x.type() == LinkType.FACING && x.hasAssembly(assembly));
                }
            }

            if(inPhaseSet && assembly.outcome() == UNSET)
                assembly.setOutcome(LINKED);

            // trim ref bases unless the assembly has been matched to a facing assembly
            if(!inFacingLink)
                RefBaseExtender.trimAssemblyRefBases(assembly, ASSEMBLY_REF_BASE_MAX_GAP);

            if(assembly.outcome() == DUP_BRANCHED)
            {
                // remove any branched assemblies which did not form a facing link
                if(inPhaseSet && inFacingLink)
                {
                    // set outcome to original assembly
                    JunctionAssembly originalAssembly = findMatchingAssembly(assembly, false);

                    if(originalAssembly != null)
                        assembly.setOutcome(originalAssembly.outcome());
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

    private boolean hasHighAssemblyCount() { return mAssemblies.size() >= HIGH_ASSEMBLY_COUNT; }

    private void checkLogPerfTime(final String stage)
    {
        if(mPerfLogTime == 0 || !hasHighAssemblyCount())
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

            SV_LOGGER.debug(format("pgId(%d) assemblies(%d: %s) stage(%s) time(%.3fs) details(links=%d candidates=%d line=%d) remoteRef(slices=%d reads=%d)",
                    mPhaseGroup.id(), mAssemblies.size(), sj, stage, seconds, mSplitLinks.size(), mExtensionCandidates.size(),
                    mLineRelatedAssemblies.size(), mRemoteRegionAssembler.remoteReadSlices(), mRemoteRegionAssembler.remoteReadsSearch()));
        }

        mStartTimeMs = System.currentTimeMillis();
    }

    // currently unused
    private void findRemoteRefCandidates()
    {
        if(!AssemblyConfig.RunRemoteRefLinking)
            return;

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
                extensionCandidate.SupportCount = supportCount + candidateCount;
                extensionCandidate.ExtraInfo = format("readSpan(%d)", remoteAssembly.refBaseLength());

                mExtensionCandidates.add(extensionCandidate);
            }
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
            applySplitLinkSupport(initialAssembly, remoteAssembly, true, true);
            initialAssembly.setOutcome(REMOTE_LINK);
            mSplitLinks.add(remoteLink);
        }
        else
        {
            applySplitLinkSupport(initialAssembly, remoteAssembly, false, true);
            mSecondarySplitLinks.add(remoteLink);
        }
    }
}
