package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_EXTENSION_BASE_MISMATCH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_REF_BASE_MAX_GAP;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_REF_SIDE_OVERLAP_BASES;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC;
import static com.hartwig.hmftools.esvee.AssemblyConstants.REF_SIDE_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.DUP_BRANCHED;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.findIndelExtensions;
import static com.hartwig.hmftools.esvee.assembly.RemoteRegionFinder.findRemoteRegions;
import static com.hartwig.hmftools.esvee.assembly.types.RefSideSoftClip.purgeRefSideSoftClips;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.CANDIDATE_DISCORDANT;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.DISCORDANT;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION_MATE;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.isDiscordantFragment;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.RefBaseAssembly;
import com.hartwig.hmftools.esvee.assembly.types.RefSideSoftClip;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;
import com.hartwig.hmftools.esvee.assembly.read.Read;

public class RefBaseExtender
{
    public RefBaseExtender() { }

    public void findAssemblyCandidateExtensions(final JunctionAssembly assembly, final List<Read> unfilteredNonJunctionReads)
    {
        // first establish potential boundaries for extending the assembly on the non-junction side
        if(assembly.indel())
        {
            // add junction mates only, could consider add reads
            findIndelExtensions(assembly, unfilteredNonJunctionReads);
            return;
        }

        int minAlignedPosition = assembly.minAlignedPosition();
        int maxAlignedPosition = assembly.maxAlignedPosition();

        boolean isForwardJunction = assembly.junction().isForward();
        int junctionPosition = assembly.junction().Position;

        // not keeping reads with unmapped mates since not sure how to incorporate their bases
        List<Read> discordantReads = unfilteredNonJunctionReads.stream()
                .filter(x -> isDiscordantCandidate(x, isForwardJunction, junctionPosition) || x.isMateUnmapped())
                .filter(x -> !assembly.hasReadSupport(x.mateRead()))
                .collect(Collectors.toList());

        List<NonJunctionRead> candidateReads = discordantReads.stream()
                .map(x -> new NonJunctionRead(x, CANDIDATE_DISCORDANT)).collect(Collectors.toList());

        discordantReads.stream().filter(x -> x.isMateUnmapped() && x.mateRead() != null).forEach(x -> assembly.addUnmappedRead(x.mateRead()));

        List<Read> remoteJunctionMates = Lists.newArrayList();
        List<Read> suppJunctionReads = Lists.newArrayList();

        // add any junction mates in the same window
        for(SupportRead support : assembly.support())
        {
            if(support.cachedRead().hasSupplementary())
                suppJunctionReads.add(support.cachedRead());

            if(isDiscordantFragment(support.cachedRead()))
            {
                remoteJunctionMates.add(support.cachedRead());
                continue;
            }

            // look to extend from local mates on the ref side of the junction
            Read mateRead = support.cachedRead().mateRead();

            if(mateRead == null || discordantReads.contains(mateRead))
                continue;

            if(!mateRead.isUnmapped())
            {
                if(isForwardJunction)
                {
                    if(mateRead.alignmentEnd() >= junctionPosition)
                        continue;
                }
                else
                {
                    if(mateRead.alignmentStart() <= junctionPosition)
                        continue;
                }

                candidateReads.add(new NonJunctionRead(mateRead, JUNCTION_MATE));
            }
            else
            {
                assembly.addUnmappedRead(mateRead);
            }
        }

        // for now add these candidate discordant reads without any further checks
        for(NonJunctionRead read : candidateReads)
        {
            assembly.addCandidateSupport(read.read(), read.type());
        }

        // now all possible discordant and junction mate reads have been collected, so test for overlaps with the min/max aligned position
        // process in order of closest to furthest-out reads in the ref base direction
        List<NonJunctionRead> sortedCandidateReads = candidateReads.stream()
                .sorted(Comparator.comparingInt(x -> isForwardJunction ? -x.read().alignmentEnd() : x.read().alignmentStart()))
                .collect(Collectors.toList());

        // look for evidence of a soft-clip on the ref side for possible branching during phasing and linking, searching from those
        // reads closest to the junction and stopping whenever a gap is encountered
        for(NonJunctionRead nonJuncRead : sortedCandidateReads)
        {
            Read read = nonJuncRead.read();

            if(isForwardJunction)
            {
                if(read.alignmentEnd() < minAlignedPosition + ASSEMBLY_REF_SIDE_OVERLAP_BASES)
                    break;

                minAlignedPosition = min(minAlignedPosition, read.alignmentStart());
            }
            else
            {
                if(read.alignmentStart() > maxAlignedPosition - ASSEMBLY_REF_SIDE_OVERLAP_BASES)
                    break;

                maxAlignedPosition = max(maxAlignedPosition, read.alignmentEnd());
            }

            assembly.checkAddRefSideSoftClip(read);
        }

        findRemoteRegions(assembly, discordantReads, remoteJunctionMates, suppJunctionReads);

        // only keep possible alternative ref-base assemblies with sufficient evidence and length
        int nonSoftClipRefPosition = isForwardJunction ? minAlignedPosition : maxAlignedPosition;
        purgeRefSideSoftClips(assembly.refSideSoftClips(), PRIMARY_ASSEMBLY_MIN_READ_SUPPORT, REF_SIDE_MIN_SOFT_CLIP_LENGTH, nonSoftClipRefPosition);
    }

    private class NonJunctionRead
    {
        private final Read mRead;
        private final SupportType mType;

        public NonJunctionRead(final Read read, final SupportType type)
        {
            mRead = read;
            mType = type;
        }

        public Read read() { return mRead; }
        public SupportType type() { return mType; }

        public String toString()
        {
            return format("%s: %s", type(), mRead);
        }
    }

    public static boolean isValidSupportCoordsVsJunction(final Read read, boolean isForwardJunction, int junctionPosition)
    {
        // cannot cross the junction since will already have considered all junction candidate reads
        // and must read in the direction of the junction
        if(isForwardJunction)
        {
            if(read.negativeStrand())
                return false;

            if(read.alignmentEnd() > junctionPosition)
                return false;
        }
        else
        {
            if(read.positiveStrand())
                return false;

            if(read.alignmentStart() < junctionPosition)
                return false;
        }

        return true;
    }

    private boolean isDiscordantCandidate(final Read read, boolean isForwardJunction, int junctionPosition)
    {
        return isValidSupportCoordsVsJunction(read, isForwardJunction, junctionPosition) && isDiscordantFragment(read);
    }

    public static void extendRefBases(
            final JunctionAssembly assembly, final List<SupportRead> candidateSupport, final RefGenomeInterface refGenome,
            boolean allowBranching)
    {
        if(candidateSupport.isEmpty())
            return;

        // find the maximal ref base extension point and make note of any recurrent soft-clip points including possible branched assemblies
        int minAlignedPosition = assembly.minAlignedPosition();
        int maxAlignedPosition = assembly.maxAlignedPosition();

        boolean isForwardJunction = assembly.junction().isForward();
        int initialRefPosition = isForwardJunction ? minAlignedPosition : maxAlignedPosition;

        // build out the reads starting with those closest to the junction, and test for gaps in ref bases
        Collections.sort(candidateSupport,
                Comparator.comparingInt(x -> isForwardJunction ? -x.cachedRead().alignmentEnd() : x.cachedRead().alignmentStart()));

        // capture RSSC from these new candidate reads
        List<RefSideSoftClip> refSideSoftClips = assembly.refSideSoftClips();
        List<SupportRead> nonJunctionSupport = Lists.newArrayListWithExpectedSize(candidateSupport.size());

        for(SupportRead support : candidateSupport)
        {
            if(isForwardJunction)
            {
                int refBaseGap = minAlignedPosition - support.alignmentEnd();

                if(refBaseGap > ASSEMBLY_REF_BASE_MAX_GAP)
                    break;

                minAlignedPosition = min(minAlignedPosition, support.alignmentStart());
            }
            else
            {
                int refBaseGap = support.alignmentStart() - maxAlignedPosition;

                if(refBaseGap > ASSEMBLY_REF_BASE_MAX_GAP)
                    break;

                maxAlignedPosition = max(maxAlignedPosition, support.alignmentEnd());
            }

            nonJunctionSupport.add(support);

            RefSideSoftClip.checkAddRefSideSoftClip(refSideSoftClips, assembly.junction(), support.cachedRead());
        }

        if(nonJunctionSupport.isEmpty())
            return;

        int nonSoftClipRefPosition = isForwardJunction ? minAlignedPosition : maxAlignedPosition;

        purgeRefSideSoftClips(refSideSoftClips, PRIMARY_ASSEMBLY_MIN_READ_SUPPORT, REF_SIDE_MIN_SOFT_CLIP_LENGTH, nonSoftClipRefPosition);

        if(refSideSoftClips.isEmpty())
        {
            // most common scenario
            extendAssemblyRefBases(assembly, nonSoftClipRefPosition, nonJunctionSupport, refGenome, false);
            return;
        }

        Collections.sort(refSideSoftClips, Comparator.comparingInt(x -> -x.readCount()));

        RefSideSoftClip candidateRefSideSoftClip = refSideSoftClips.get(0);

        int nonSoftClipSupport = 0;

        if(isForwardJunction)
        {
            nonSoftClipSupport = (int)nonJunctionSupport.stream()
                    .filter(x -> !x.isLeftClipped()).filter(x -> x.alignmentStart() < candidateRefSideSoftClip.Position).count();
        }
        else
        {
            nonSoftClipSupport = (int)nonJunctionSupport.stream()
                    .filter(x -> !x.isRightClipped()).filter(x -> x.alignmentEnd() > candidateRefSideSoftClip.Position).count();
        }

        // select the ref position with the most support
        int primaryRefPosition, primaryRefPositionSupport, secondRefPositionSupport;
        boolean usesSoftClippedPosition = false;

        if(nonSoftClipSupport > candidateRefSideSoftClip.readCount())
        {
            primaryRefPosition = nonSoftClipRefPosition;
            primaryRefPositionSupport = nonSoftClipSupport;
            secondRefPositionSupport = candidateRefSideSoftClip.readCount();
        }
        else
        {
            primaryRefPosition = candidateRefSideSoftClip.Position;
            primaryRefPositionSupport = candidateRefSideSoftClip.readCount();
            secondRefPositionSupport = nonSoftClipSupport;
            usesSoftClippedPosition = true;
        }

        double secondRefPositionSupportPerc = secondRefPositionSupport / (double)primaryRefPositionSupport;
        boolean hasSufficientSecondRefSupport = secondRefPositionSupport >= PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT
                && secondRefPositionSupportPerc >= PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC;

        if(!allowBranching || !hasSufficientSecondRefSupport)
        {
            extendAssemblyRefBases(assembly, primaryRefPosition, nonJunctionSupport, refGenome, usesSoftClippedPosition);
            return;
        }

        branchAssembliesFromRefBases(assembly, nonSoftClipRefPosition, nonJunctionSupport, primaryRefPositionSupport, refGenome);
    }

    private static void extendAssemblyRefBases(
            final JunctionAssembly assembly, int newRefPosition, final List<SupportRead> nonJunctionSupport,
            final RefGenomeInterface refGenome, boolean isSoftClipped)
    {
        int initialRefPosition = assembly.isForwardJunction() ? assembly.minAlignedPosition() : assembly.maxAlignedPosition();
        boolean hasPriorRefExtension = assembly.support().stream().anyMatch(x -> x.type() == JUNCTION_MATE || x.type() == DISCORDANT);

        if(newRefPosition == initialRefPosition || hasPriorRefExtension)
        {
            // add in supporting reads
            for(SupportRead support : nonJunctionSupport)
            {
                assembly.addDiscordantSupport(support, ASSEMBLY_EXTENSION_BASE_MISMATCH);
                support.clearCachedRead();
            }

            return;
        }

        if(assembly.isForwardJunction() == (newRefPosition < initialRefPosition))
        {
            RefBaseAssembly refBaseAssembly = new RefBaseAssembly(assembly, newRefPosition, refGenome);

            checkAddRefAssemblySupport(refBaseAssembly, nonJunctionSupport, Collections.emptySet());

            // SV_LOGGER.debug("assembly({}) post-support bases: {}", refBaseAssembly, new String(refBaseAssembly.bases()));

            if(refBaseAssembly.supportCount() > 0)
                assembly.mergeRefBaseAssembly(refBaseAssembly, assembly.refBaseLength(), "non-branched");
        }
        else if(isSoftClipped)
        {
            // need to trim the existing ref bases due to ref-side soft-clipped determining the bounds
            assembly.trimRefBases(newRefPosition);
        }
    }

    private static void branchAssembliesFromRefBases(
            final JunctionAssembly originalAssembly, int nonSoftClipRefPosition,
            final List<SupportRead> nonJunctionSupport, int maxRefSideSupport, final RefGenomeInterface refGenome)
    {
        // now build out any distinct, branching assemblies from ref-base soft-clips
        List<Set<String>> excludedReadIdsList = allocateExcludedReads(originalAssembly, nonJunctionSupport);
        List<SupportRead> initialSupport = Lists.newArrayList(originalAssembly.support());

        List<JunctionAssembly> branchedAssemblies = Lists.newArrayList(originalAssembly);

        // the original assembly is handled first
        for(int i = 0; i < excludedReadIdsList.size(); ++i)
        {
            Set<String> excludedReads = excludedReadIdsList.get(i);

            JunctionAssembly junctionAssembly = null;

            int extensionRefPosition;

            if(i == 0)
            {
                junctionAssembly = originalAssembly;

                // first remove support from the main assembly
                originalAssembly.removeSupportReads(excludedReads);
                extensionRefPosition = nonSoftClipRefPosition;
            }
            else
            {
                RefSideSoftClip refSideSoftClip = originalAssembly.refSideSoftClips().get(i - 1);

                if(refSideSoftClip.matchesOriginal() || refSideSoftClip.hasProximateMatch(originalAssembly.refBasePosition()))
                    continue;

                int refBaseDifference = refSideSoftClip.Position - originalAssembly.refBasePosition();

                if(originalAssembly.junction().isReverse())
                    refBaseDifference *= -1;

                if(refBaseDifference <= 0)
                    break;

                int newRefBaseLength = originalAssembly.refBaseLength() - refBaseDifference;

                junctionAssembly = new JunctionAssembly(originalAssembly, refSideSoftClip, newRefBaseLength, initialSupport, excludedReads);

                extensionRefPosition = refSideSoftClip.Position;
            }

            RefBaseAssembly refBaseAssembly = new RefBaseAssembly(junctionAssembly, extensionRefPosition, refGenome);

            checkAddRefAssemblySupport(refBaseAssembly, nonJunctionSupport, excludedReads);

            // only add branched assemblies if they have sufficient support
            if(junctionAssembly == originalAssembly)
            {
                if(refBaseAssembly.supportCount() > 0) // same criteria as above
                {
                    junctionAssembly.mergeRefBaseAssembly(refBaseAssembly, junctionAssembly.refBaseLength(),"branched original");
                }
            }
            else
            {
                if(refBaseAssembly.supportCount() == 0)
                    continue;

                junctionAssembly.mergeRefBaseAssembly(refBaseAssembly, originalAssembly.refBaseLength(), "branched new");

                if(AssemblyUtils.hasUnsetBases(junctionAssembly))
                    continue;

                // same criteria as above
                double supportPerc = junctionAssembly.supportCount() / (double)maxRefSideSupport;
                boolean hasSufficientSecondRefSupport = junctionAssembly.supportCount() >= PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT
                        && supportPerc >= PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC;

                if(!hasSufficientSecondRefSupport)
                    continue;

                junctionAssembly.buildRepeatInfo();
                branchedAssemblies.add(junctionAssembly);
                junctionAssembly.setOutcome(DUP_BRANCHED);

                // for now only 1 branched assembly will be made
                break;
            }
        }

        // set references between them - for now just for TSV output
        for(JunctionAssembly junctionAssembly : branchedAssemblies)
        {
            if(junctionAssembly != originalAssembly)
                originalAssembly.phaseGroup().addDerivedAssembly(junctionAssembly);
        }
    }

    private static void checkAddRefAssemblySupport(
            final RefBaseAssembly refBaseAssembly, final List<SupportRead> nonJunctionReads, final Set<String> excludedReadIds)
    {
        List<SupportRead> newSupportReads = nonJunctionReads.stream().filter(x -> !excludedReadIds.contains(x.id())).collect(Collectors.toList());

        // favour junction mates first
        Collections.sort(newSupportReads, Comparator.comparingInt(x -> x.type() == JUNCTION_MATE ? 0 : 1));

        for(SupportRead support : newSupportReads)
        {
            SupportType type = support.type() == CANDIDATE_DISCORDANT ? DISCORDANT : support.type();
            refBaseAssembly.checkAddRead(support.cachedRead(), type, ASSEMBLY_EXTENSION_BASE_MISMATCH, ASSEMBLY_REF_SIDE_OVERLAP_BASES);
        }
    }

    private static List<Set<String>> allocateExcludedReads(final JunctionAssembly assembly, final List<SupportRead> nonJunctionReads)
    {
        // now build out any distinct, branching assemblies from ref-base soft-clips
        List<RefSideSoftClip> refSideSoftClips = assembly.refSideSoftClips();

        boolean isForwardJunction = assembly.isForwardJunction();
        Set<String> nscReads = Sets.newHashSet();

        int maxRefSideSoftClipPosition = isForwardJunction ?
                refSideSoftClips.stream().filter(x -> !x.matchesOriginal()).mapToInt(x -> x.Position).min().orElse(0)
                : refSideSoftClips.stream().filter(x -> !x.matchesOriginal()).mapToInt(x -> x.Position).max().orElse(0);

        List<SupportRead> allReads = Lists.newArrayList(assembly.support());
        allReads.addAll(nonJunctionReads);

        Set<String> softClippedReadIds = Sets.newHashSet();
        refSideSoftClips.stream().filter(x -> !x.matchesOriginal()).forEach(x -> softClippedReadIds.addAll(x.readIds()));

        // find any supporting read or read mate which extends past the furthest soft-clip position
        for(SupportRead read : allReads)
        {
            if(softClippedReadIds.stream().anyMatch(x -> x.equals(read.id())))
                continue;

            boolean addRead = false;

            if(isForwardJunction)
            {
                if(read.alignmentStart() < maxRefSideSoftClipPosition
                || (!read.isDiscordant() && read.isMateMapped() && read.mateAlignmentStart() < maxRefSideSoftClipPosition))
                {
                    addRead = true;
                }
            }
            else
            {
                if(read.alignmentEnd() > maxRefSideSoftClipPosition
                || (!read.isDiscordant() && read.isMateMapped() && read.mateAlignmentEnd() > maxRefSideSoftClipPosition))
                {
                    addRead = true;
                }
            }

            if(addRead)
            {
                nscReads.add(read.id());
            }
        }

        int finalAssemblyCount = 1 + refSideSoftClips.size();
        List<Set<String>> excludedReadsList = Lists.newArrayListWithCapacity(finalAssemblyCount);

        // the original assembly cannot have any
        excludedReadsList.add(softClippedReadIds);

        for(RefSideSoftClip refSideSoftClip : refSideSoftClips)
        {
            Set<String> excludedReads = Sets.newHashSet();

            if(!refSideSoftClip.matchesOriginal())
                excludedReads.addAll(nscReads);

            excludedReadsList.add(excludedReads);

            for(RefSideSoftClip other : refSideSoftClips)
            {
                if(refSideSoftClip == other)
                    continue;

                excludedReads.addAll(other.readIds());
            }
        }

        return excludedReadsList;
    }
}
