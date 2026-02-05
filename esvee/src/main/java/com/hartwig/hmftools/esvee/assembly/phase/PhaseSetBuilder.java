package com.hartwig.hmftools.esvee.assembly.phase;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_REF_BASE_MAX_GAP;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.LOCAL_ASSEMBLY_MATCH_DISTANCE;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PROXIMATE_DUP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.RefBaseExtender.checkAddRefBaseRead;
import static com.hartwig.hmftools.esvee.assembly.phase.AssemblyLinker.isAssemblyIndelLink;
import static com.hartwig.hmftools.esvee.assembly.phase.AssemblyLinker.isFacingAssemblyCandidate;
import static com.hartwig.hmftools.esvee.assembly.phase.ExtensionType.INDEL;
import static com.hartwig.hmftools.esvee.assembly.phase.ExtensionType.LOCAL_DEL_DUP;
import static com.hartwig.hmftools.esvee.assembly.phase.ExtensionType.SPLIT_LINK;
import static com.hartwig.hmftools.esvee.assembly.phase.ExtensionType.UNMAPPED;
import static com.hartwig.hmftools.esvee.assembly.phase.RemoteReadExtractor.collectCandidateRemoteRegions;
import static com.hartwig.hmftools.esvee.assembly.phase.RemoteReadExtractor.purgeSupplementaryReads;
import static com.hartwig.hmftools.esvee.assembly.read.Read.findMatchingFragmentSupport;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.isDiscordantFragment;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyLink.swapAssemblies;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.DUP_BRANCHED;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.LINKED;
import static com.hartwig.hmftools.esvee.assembly.RefBaseExtender.extendRefBases;
import static com.hartwig.hmftools.esvee.assembly.phase.AssemblyLinker.tryAssemblyFacing;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.LOCAL_INDEL;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.SECONDARY;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.UNSET;
import static com.hartwig.hmftools.esvee.assembly.types.SupportRead.hasFragmentOtherRead;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.DISCORDANT;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.EXTENSION;
import static com.hartwig.hmftools.esvee.common.CommonUtils.isLineInsertPair;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.hasPairedReads;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.sv.SvUtils;
import com.hartwig.hmftools.common.utils.Integers;
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

import org.jetbrains.annotations.Nullable;

public class PhaseSetBuilder
{
    private final PhaseGroup mPhaseGroup;
    private final RefGenomeInterface mRefGenome;
    private final RemoteReadExtractor mRemoteReadExtractor;
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
    private final boolean mHasHighAssemblyCount;
    private final boolean mHasUltraHighAssemblyCount;
    private long mStartTimeMs;
    private boolean mDetailsLogged;

    private Stage mCurrentStage;
    private int mRoutineIteration;

    private enum Stage
    {
        FormLocalLinks,
        FindLineExtensions,
        FindUnmappedExtensions,
        FindNonLocalLinkCandidates,
        FindNonLocalLinkCandidatesPostExtensions,
        ApplyLinksAndExtensions,
        AddUnlinkedAssemblyRefSupport,
        FinalisePhaseSets,
        MergePhaseSetAlignments;
    }

    public PhaseSetBuilder(
            final RefGenomeInterface refGenome, final RemoteReadExtractor remoteReadExtractor, final PhaseGroup phaseGroup)
    {
        mRefGenome = refGenome;
        mPhaseGroup = phaseGroup;

        mRemoteReadExtractor = remoteReadExtractor;
        mRemoteReadExtractor.resetCounts();

        int assemblyCount = mPhaseGroup.assemblies().size();
        mHasHighAssemblyCount = assemblyCount >= HIGH_ASSEMBLY_COUNT;
        mHasUltraHighAssemblyCount = assemblyCount >= ULTRA_HIGH_ASSEMBLY_COUNT;

        if(mHasUltraHighAssemblyCount)
        {
            mAssemblies = filterUltraHighAssemblies(mPhaseGroup);
        }
        else
        {
            mAssemblies = mPhaseGroup.assemblies();
        }

        mLocalSequenceMatcher = new LocalSequenceMatcher(refGenome, LOCAL_ASSEMBLY_MATCH_DISTANCE);

        mPhaseSets = mPhaseGroup.phaseSets();
        mSecondarySplitLinks = mPhaseGroup.secondaryLinks();

        mLineRelatedAssemblies = Sets.newHashSet();

        mExtensionCandidates = Lists.newArrayList();
        mSplitLinks = Lists.newArrayList();
        mFacingLinks = Lists.newArrayList();
        mBranchedAssemblies = Lists.newArrayList();
        mLocallyLinkedAssemblies = Sets.newHashSet();

        mStartTimeMs = 0;
        mRoutineIteration = 0;
        mCurrentStage = null;
        mDetailsLogged = false;
    }

    private static final int HIGH_ASSEMBLY_COUNT = 100;
    private static final int ULTRA_HIGH_ASSEMBLY_COUNT = 1000;
    private static final int HIGH_ASSEMBLY_READ_COUNT = 100;
    private static final int MAX_HIGH_ASSEMBLY_MATCHED_READS = 1000;

    public void run()
    {
        buildPhaseSets();

        initialisePerfStage(Stage.MergePhaseSetAlignments);

        mPhaseGroup.finalisePhaseSetAlignments();
        checkLogPerfTime();
    }

    public void buildPhaseSets()
    {
        findLineExtensions();

        findLocalLinks();

        findOtherLinksAndExtensions();

        addUnlinkedAssemblyRefSupport();

        initialisePerfStage(Stage.FinalisePhaseSets);

        formFacingLinks();

        checkBranchedAssemblies();

        formPhaseSets();

        addChainedSupport();

        cleanupAssemblies();

        checkLogPerfTime(); // covers facing links onwards
    }

    private void findLocalLinks()
    {
        initialisePerfStage(Stage.FormLocalLinks);

        // find local candidate links
        findSplitLinkCandidates(true, null);

        // prioritise and capture local links
        Collections.sort(mExtensionCandidates, new ExtensionCandidate.LocalLinkComparator());

        for(ExtensionCandidate extensionCandidate : mExtensionCandidates)
        {
            if(!extensionCandidate.isValid())
                continue;

            AssemblyLink assemblyLink = extensionCandidate.Link;

            extensionCandidate.markSelected();

            boolean isPrimaryLink = !mLocallyLinkedAssemblies.contains(extensionCandidate.Assembly)
                    && !mLocallyLinkedAssemblies.contains(extensionCandidate.SecondAssembly);

            applySplitLink(assemblyLink, isPrimaryLink, true);

            mLocallyLinkedAssemblies.add(extensionCandidate.Assembly);
            mLocallyLinkedAssemblies.add(extensionCandidate.SecondAssembly);
        }

        // check for other local alignments

        for(JunctionAssembly assembly : mAssemblies)
        {
            if(!mLocallyLinkedAssemblies.contains(assembly))
            {
                formsLocalLink(assembly);
            }
        }

        checkLogPerfTime();

        // no need to keep this info for subsequent non-local candidate analysis
        mExtensionCandidates.clear();
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

    private void findSplitLinkCandidates(boolean localIndelOnly, @Nullable final Set<JunctionAssembly> extendedAssemblies)
    {
        if(mAssemblies.size() < 2)
            return;

        // test each assembly pair - restricted to local-only links if specified
        // if not, check that a pair hasn't already been tested
        // if no support is found then no need to check the link

        for(int i = 0; i < mAssemblies.size() - 1; ++i)
        {
            JunctionAssembly assembly1 = mAssemblies.get(i);

            if(!localIndelOnly && !allowLocalIndelSecondaryLinks(assembly1)) // no attempt to find secondaries for these
                continue;

            Set<String> firstReadIds = Sets.newHashSet();
            boolean firstHighReadCount = assemblyHasHighReadCount(assembly1);

            if(!localIndelOnly)
                populateReadIds(assembly1, firstReadIds, localIndelOnly);

            // allow linking assemblies to be included, to allow secondary links to be found
            for(int j = i + 1; j < mAssemblies.size(); ++j)
            {
                JunctionAssembly assembly2 = mAssemblies.get(j);

                if(!localIndelOnly && !allowLocalIndelSecondaryLinks(assembly2))
                    continue;

                if(extendedAssemblies != null)
                {
                    // ignore if neither assembly has been extended
                    if(!extendedAssemblies.contains(assembly1) && !extendedAssemblies.contains(assembly2))
                        continue;

                    // or if the link has already been formed
                    if(mSplitLinks.stream().anyMatch(x -> x.matches(assembly1, assembly2)))
                        continue;
                }

                boolean isLocalIndel = isAssemblyIndelLink(assembly1, assembly2);
                boolean isLocalLinkCandidate = false;

                if(!assembly1.discordantOnly() && !assembly2.discordantOnly() && assembly1.indel() == assembly2.indel())
                {
                    isLocalLinkCandidate = isLocalIndel || isLocalAssemblyCandidate(assembly1, assembly2);
                }

                if(localIndelOnly)
                {
                    if(!isLocalLinkCandidate)
                        continue;
                }
                else
                {
                    if(isLocalLinkCandidate) // avoid checking the same pair again, since they were checked in the local-only call
                        continue;
                }

                Set<String> secondReadIds = Sets.newHashSet();
                boolean secondHighReadCount = assemblyHasHighReadCount(assembly2);

                if(!localIndelOnly)
                    populateReadIds(assembly2, secondReadIds, localIndelOnly);

                Set<String> firstReadIdsRelated;
                Set<String> secondReadIdsRelated;

                if(!localIndelOnly && (firstHighReadCount || secondHighReadCount))
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
                boolean hasSharedFragments = localIndelOnly || hasSharedFragment(firstReadIdsRelated, secondReadIdsRelated);

                AssemblyLink assemblyLink = null;

                if(hasSharedFragments)
                    assemblyLink = checkSplitLink(assembly1, assembly2, isLocalLinkCandidate);

                if(!hasSharedFragments || assemblyLink == null)
                {
                    if(localIndelOnly)
                    {
                        // cache to avoid checking on a second pass
                        mExtensionCandidates.add(new ExtensionCandidate(LOCAL_DEL_DUP, assembly1, assembly2));
                    }

                    continue;
                }

                ExtensionType type = isLocalLinkCandidate ? (isLocalIndel ? INDEL : LOCAL_DEL_DUP) : SPLIT_LINK;

                ExtensionCandidate extensionCandidate = new ExtensionCandidate(type, assemblyLink);
                mExtensionCandidates.add(extensionCandidate);

                if(localIndelOnly)
                {
                    // populate read IDs now that the link has been made
                    if(firstReadIds.isEmpty())
                        populateReadIds(assembly1, firstReadIds, localIndelOnly);

                    populateReadIds(assembly2, secondReadIds, localIndelOnly);
                }

                // now count up all possible linking fragments to compare with other candidate links and extensions
                countSharedFragments(extensionCandidate, firstReadIdsRelated, secondReadIdsRelated);
            }

            if(checkMaxRoutineIteration())
                return;
        }
    }

    private boolean allowLocalIndelSecondaryLinks(final JunctionAssembly assembly)
    {
        // protect local indels from testing and making secondary links, with exceptions made for those at LINE insertion sites
        // or others with long insertion sequences
        if(mLineRelatedAssemblies.contains(assembly))
            return true;

        if(!mLocallyLinkedAssemblies.contains(assembly))
            return true;

        AssemblyLink existingLink = mSplitLinks.stream().filter(x -> x.hasAssembly(assembly)).findFirst().orElse(null);

        return existingLink == null || existingLink.insertedBases().length() > existingLink.length();
    }

    private static boolean assemblyHasHighReadCount(final JunctionAssembly assembly)
    {
        return assembly.supportCount() >= HIGH_ASSEMBLY_READ_COUNT || assembly.candidateSupport().size() >= HIGH_ASSEMBLY_READ_COUNT;
    }

    private void populateReadIds(final JunctionAssembly assembly, final Set<String> readIds, boolean localIndelOnly)
    {
        if(assemblyHasHighReadCount(assembly))
        {
            if(localIndelOnly)
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

            // ignore candidate discordant reads around local indels since it's not clear enough if it relates to the assembly
            if(!localIndelOnly)
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
        // check non-local links after first attempting to extend bases using unmapped and remote reads
        findUnmappedExtensions();

        Set<JunctionAssembly> extendedAssemblies = applyOtherLinksAndExtensions(Stage.FindNonLocalLinkCandidates, null);

        if(!extendedAssemblies.isEmpty())
        {
            // if some assemblies were extended, they may now form links when they didn't before - so check this subset again
            applyOtherLinksAndExtensions(Stage.FindNonLocalLinkCandidatesPostExtensions, extendedAssemblies);
        }
    }

    private Set<JunctionAssembly> applyOtherLinksAndExtensions(final Stage stage, @Nullable final Set<JunctionAssembly> extendedAssemblies)
    {
        initialisePerfStage(stage);

        findSplitLinkCandidates(false, extendedAssemblies); // since local candidate links have already been found and applied

        checkLogPerfTime();

        initialisePerfStage(Stage.ApplyLinksAndExtensions);

        // on the second pass only check asemblies which have been extended for new split links
        boolean checkOnlyExtendedAssemblies = extendedAssemblies != null;

        // prioritise and select from all remaining candidates
        List<ExtensionCandidate> remainingCandidates = mExtensionCandidates.stream()
                .filter(x -> !x.selected()) // skip those already registered
                .filter(x -> x.isValid())
                .filter(x -> !checkOnlyExtendedAssemblies || x.Type != UNMAPPED)
                .collect(Collectors.toList());

        if(remainingCandidates.isEmpty())
            return Collections.emptySet();

        Collections.sort(remainingCandidates, new ExtensionCandidate.StandardComparator());

        Set<JunctionAssembly> primaryLinkedAssemblies = Sets.newHashSet(mLocallyLinkedAssemblies);
        mSplitLinks.forEach(x -> primaryLinkedAssemblies.add(x.first()));
        mSplitLinks.forEach(x -> primaryLinkedAssemblies.add(x.second()));

        Set<JunctionAssembly> newlyExtendedAssemblies = Sets.newHashSet();

        for(ExtensionCandidate extensionCandidate : remainingCandidates)
        {
            if(extensionCandidate.Type == SPLIT_LINK)
            {
                if(checkOnlyExtendedAssemblies)
                {
                    if(!extendedAssemblies.contains(extensionCandidate.Assembly)
                    && !extendedAssemblies.contains(extensionCandidate.SecondAssembly))
                    {
                        continue;
                    }
                }

                boolean eitherInPrimary = primaryLinkedAssemblies.contains(extensionCandidate.Assembly)
                        || primaryLinkedAssemblies.contains(extensionCandidate.SecondAssembly);

                if(eitherInPrimary
                && (extensionCandidate.Assembly.discordantOnly() || extensionCandidate.SecondAssembly.discordantOnly()))
                {
                    // ignore discordant-only junctions with links to an assembly already link
                    continue;
                }

                extensionCandidate.markSelected();
                applySplitLink(extensionCandidate.Link, !eitherInPrimary, false);

                if(!eitherInPrimary)
                {
                    primaryLinkedAssemblies.add(extensionCandidate.Assembly);
                    primaryLinkedAssemblies.add(extensionCandidate.SecondAssembly);
                }
            }
            else if(extensionCandidate.Type == UNMAPPED)
            {
                // extend the assembly if not already linked in any way - this facilitates further links and/or a better SGL for alignment
                if(extensionCandidate.Assembly.outcome() == UNSET)
                {
                    extensionCandidate.markSelected();
                    applyUnmappedReadExtension(extensionCandidate);
                    newlyExtendedAssemblies.add(extensionCandidate.Assembly);
                }
            }
        }

        checkLogPerfTime();

        return newlyExtendedAssemblies;
    }

    private void findUnmappedExtensions()
    {
        if(!hasPairedReads())
            return;

        initialisePerfStage(Stage.FindUnmappedExtensions);

        double extractRemoteReadsTotalSeconds = 0;

        int totalRawRemoteRegions = 0;
        int totalRemoteRegions = 0;

        // any assembly not in a link uses unmapped reads to try to extend the extension sequence
        for(JunctionAssembly assembly : mAssemblies)
        {
            if(mLocallyLinkedAssemblies.contains(assembly) || mLineRelatedAssemblies.contains(assembly)) // ignore if already processed as a line site or a local indel
                continue;

            List<Read> unmappedReads = Lists.newArrayList(assembly.unmappedReads());

            long startTimeMs = System.currentTimeMillis();
            List<RemoteRegion> remoteRegions = collectCandidateRemoteRegions(assembly, mAssemblies, mHasHighAssemblyCount);

            totalRawRemoteRegions += assembly.remoteRegions().size();

            List<Read> remoteReads = mRemoteReadExtractor.extractRemoteRegionReads(mPhaseGroup.id(), remoteRegions, mHasHighAssemblyCount);
            extractRemoteReadsTotalSeconds += (System.currentTimeMillis() - startTimeMs) / 1000.0;

            totalRemoteRegions += remoteRegions.size();

            purgeSupplementaryReads(assembly, remoteReads);
            unmappedReads.addAll(remoteReads);

            if(unmappedReads.isEmpty())
                break;

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

            if(checkMaxRoutineIteration())
                break;
        }

        checkLogPerfTime();

        if(AssemblyConfig.PerfLogTime > 0 && extractRemoteReadsTotalSeconds >= AssemblyConfig.PerfLogTime)
        {
            SV_LOGGER.debug(format("%s phase set stage(%s) remoteRegions(raw=%d collected=%d slices=%d) reads(%d) extractReads(%.1f)",
                    getPhaseGroupInfo(), mCurrentStage, totalRawRemoteRegions, totalRemoteRegions, mRemoteReadExtractor.remoteReadSlices(),
                    mRemoteReadExtractor.remoteReadsSearch(), extractRemoteReadsTotalSeconds));
        }
    }

    private void findLineExtensions()
    {
        // find line insertion sites - pairs of proximate INDEL-type junctions with a line motif - and jointly search
        // for extension reads across the two assemblies
        mAssemblies.stream().filter(x -> x.hasLineSequence()).forEach(x -> mLineRelatedAssemblies.add(x));

        if(mLineRelatedAssemblies.isEmpty())
            return;

        initialisePerfStage(Stage.FindLineExtensions);

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

            if(!hasPairedReads())
                continue;

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

            List<Read> remoteReads = mRemoteReadExtractor.extractRemoteRegionReads(mPhaseGroup.id(), combinedRemoteRegions, mHasHighAssemblyCount);
            purgeSupplementaryReads(assembly, remoteReads);

            sharedUnmappedReads.addAll(remoteReads);

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

        checkLogPerfTime();
    }

    private void applySplitLink(final AssemblyLink assemblyLink, boolean isPrimaryLink, boolean isLocalIndel)
    {
        boolean allowBranching = isPrimaryLink && !(assemblyLink.svType() == DUP && assemblyLink.length() < PROXIMATE_DUP_LENGTH);

        // discordant reads are not used to build out local indel links
        boolean allowDiscordantReads = !isLocalIndel && !SvUtils.isShortLocalDelDupIns(assemblyLink.svType(), assemblyLink.length());

        if(isLocalIndel)
        {
            adjustLocalIndelSupport(assemblyLink.first(), assemblyLink.second());
        }

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
        initialisePerfStage(Stage.AddUnlinkedAssemblyRefSupport);

        // any assembly which did not form a link or only an unmapped extension will now extend its ref bases from junction & extension mates
        List<JunctionAssembly> assemblies = Lists.newArrayList(mAssemblies);

        boolean allowRefSideSoftClipBranching = true;

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

            extendRefBases(assembly, refExtensionReads, mRefGenome, allowRefSideSoftClipBranching);
        }

        checkLogPerfTime();
    }

    private boolean applySplitLinkSupport(
            final JunctionAssembly assembly1, final JunctionAssembly assembly2, boolean allowBranching, boolean allowDiscordantReads)
    {
        if(!hasPairedReads())
            return false;

        // look for shared reads between the assemblies, and factor in discordant reads which were only considered candidates until now
        List<Read> matchedCandidates1 = Lists.newArrayList();
        List<Read> matchedCandidates2 = Lists.newArrayList();

        boolean hasHighAssemblyCount = mHasHighAssemblyCount;
        boolean candidateSupportLimited1 = hasHighAssemblyCount && assembly1.candidateSupport().size() > MAX_HIGH_ASSEMBLY_MATCHED_READS;
        boolean candidateSupportLimited2 = hasHighAssemblyCount && assembly2.candidateSupport().size() > MAX_HIGH_ASSEMBLY_MATCHED_READS;

        List<Read> candidateSupport1 = candidateSupportLimited1 ?
                assembly1.candidateSupport().subList(0, MAX_HIGH_ASSEMBLY_MATCHED_READS) : assembly1.candidateSupport();

        List<Read> candidateSupport2 = candidateSupportLimited2
                ? assembly2.candidateSupport().subList(0, MAX_HIGH_ASSEMBLY_MATCHED_READS) : assembly2.candidateSupport();

        checkMatchingCandidateSupport(
                assembly2, allowDiscordantReads, candidateSupport1, candidateSupport2, matchedCandidates1, matchedCandidates2);

        checkMatchingCandidateSupport(
                assembly1, allowDiscordantReads, candidateSupport2, Collections.emptyList(), matchedCandidates2, matchedCandidates1);

        // remove candidates from actual assembly candidate lists, since only subset lists were manipulated
        if(candidateSupportLimited1)
            matchedCandidates1.forEach(x -> assembly1.candidateSupport().remove(x));

        if(candidateSupportLimited2)
            matchedCandidates2.forEach(x -> assembly2.candidateSupport().remove(x));

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
        extendRefBases(assembly1, matchedCandidates1, mRefGenome, allowBranching);
        extendRefBases(assembly2, matchedCandidates2, mRefGenome, allowBranching);

        // register any newly branched assemblies
        for(JunctionAssembly assembly : mAssemblies)
        {
            if(assembly.outcome() == DUP_BRANCHED && !mBranchedAssemblies.contains(assembly))
                mBranchedAssemblies.add(assembly);
        }

        return true;
    }

    private void checkMatchingCandidateSupport(
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

    private void adjustLocalIndelSupport(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
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

        List<AssemblyLink> splitLinks = Lists.newArrayList(mSplitLinks);
        splitLinks.addAll(mSecondarySplitLinks);

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

                AssemblyLink facingLink = tryAssemblyFacing(assembly1, assembly2, splitLinks);

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

    private void initialisePerfStage(final Stage stage)
    {
        mCurrentStage = stage;
        mStartTimeMs = System.currentTimeMillis();
        mRoutineIteration = 0;
    }

    private void checkLogPerfTime()
    {
        if(AssemblyConfig.PerfLogTime == 0 || !mHasHighAssemblyCount)
            return;

        long timeTakenMs = System.currentTimeMillis() - mStartTimeMs;
        double seconds = timeTakenMs / 1000.0;

        if(seconds >= AssemblyConfig.PerfLogTime)
        {
            if(mCurrentStage == Stage.MergePhaseSetAlignments)
            {
                SV_LOGGER.debug(format("%s phase set stage(%s) time(%.1fs) phaseSets(%d)",
                        getPhaseGroupInfo(), mCurrentStage, seconds, mPhaseGroup.phaseSets().size()));
            }
            else
            {
                SV_LOGGER.debug(format("%s phase set stage(%s) time(%.1fs) details(links=%d candidates=%d line=%d)",
                        getPhaseGroupInfo(), mCurrentStage, seconds, mSplitLinks.size(), mExtensionCandidates.size(),
                        mLineRelatedAssemblies.size()));
            }
        }
    }

    private boolean checkMaxRoutineIteration()
    {
        if(AssemblyConfig.PhaseProcessingLimit == 0)
            return false;

        ++mRoutineIteration;

        if(mRoutineIteration < AssemblyConfig.PhaseProcessingLimit)
            return false;

        SV_LOGGER.debug(format("%s exiting phase set stage(%s) details(links=%d candidates=%d) after %d iterations",
                getPhaseGroupInfo(), mCurrentStage, mSplitLinks.size(), mExtensionCandidates.size(), mRoutineIteration));

        return true;
    }

    private String getPhaseGroupInfo()
    {
        if(mDetailsLogged)
            return format("pgId(%d)", mPhaseGroup.id());

        mDetailsLogged = true;

        StringJoiner sj = new StringJoiner(";");
        for(int i = 0; i < min(mAssemblies.size(), 4); ++i)
        {
            sj.add(mAssemblies.get(i).junction().coords());
        }

        return format("pgId(%d) assemblies(%d: %s)", mPhaseGroup.id(), mAssemblies.size(), sj);
    }

    private List<JunctionAssembly> filterUltraHighAssemblies(final PhaseGroup phaseGroup)
    {
        // log stats to make further analysis easier
        int lineCount = 0;
        int splitCount = 0;
        int indelCount = 0;
        int discordantCount = 0;

        List<Integer> extensionLengths = Lists.newArrayList();
        List<Integer> fragmentCounts = Lists.newArrayList();
        List<Integer> mismatchReadCounts = Lists.newArrayList();

        for(JunctionAssembly assembly : phaseGroup.assemblies())
        {
            if(assembly.junction().DiscordantOnly)
            {
                ++discordantCount;
            }
            else if(assembly.junction().indelBased())
            {
                ++indelCount;
            }
            else
            {
                ++splitCount;

                if(assembly.hasLineSequence())
                    ++lineCount;
            }

            fragmentCounts.add(assembly.supportCount());
            extensionLengths.add(assembly.extensionLength());
            mismatchReadCounts.add(assembly.mismatchReadCount());
        }

        int medianFragCount = (int)round(Integers.median(fragmentCounts));
        int medianExtensionLength = (int)round(Integers.median(extensionLengths));
        int medianMismatchReads = (int)round(Integers.median(mismatchReadCounts));

        SV_LOGGER.info("pgId({}) assemblies({}) types(disc={} indel={} split={} line={}) median(frags={} extLength={} mismatchReads={})",
                mPhaseGroup.id(), mPhaseGroup.assemblies().size(), discordantCount, indelCount, splitCount, lineCount,
                medianFragCount, medianExtensionLength, medianMismatchReads);

        // take the top X assemblies only, scored by longest extensions, lowest mismatch ratio and highest frag counts
        List<JunctionAssembly> allAssemblies = Lists.newArrayList(phaseGroup.assemblies());

        Collections.sort(allAssemblies, new JunctionAssemblyScoredSorter());

        return allAssemblies.subList(0, ULTRA_HIGH_ASSEMBLY_COUNT);
    }

    private static double scoreJunctionAssembly(final JunctionAssembly assembly)
    {
        int fragCount = assembly.supportCount();
        int mismatchReadCount = assembly.mismatchReadCount();
        int extBaseLength = assembly.extensionLength();
        int refBaseLength = assembly.refBaseLength();

        double fragScore = fragCount / (double)max(mismatchReadCount, 1);
        double lengthScore = extBaseLength * refBaseLength / (double)(MIN_VARIANT_LENGTH * MIN_VARIANT_LENGTH);
        return fragScore * lengthScore;
    }

    private static class JunctionAssemblyScoredSorter implements Comparator<JunctionAssembly>
    {
        public int compare(final JunctionAssembly first, final JunctionAssembly second)
        {
            double firstScore = scoreJunctionAssembly(first);
            double secondScore = scoreJunctionAssembly(second);

            return Double.compare(-firstScore, -secondScore);
        }
    }
}
