package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_EXTENSION_BASE_MISMATCH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_REF_SIDE_OVERLAP_BASES;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.AssemblyConstants.REF_SIDE_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.findUnsetBases;
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
                .filter(x -> isDiscordantCandidate(x, isForwardJunction, junctionPosition))
                .filter(x -> !assembly.hasReadSupport(x.mateRead()))
                .collect(Collectors.toList());

        List<NonJunctionRead> candidateReads = discordantReads.stream()
                .map(x -> new NonJunctionRead(x, CANDIDATE_DISCORDANT)).collect(Collectors.toList());

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

            if(mateRead == null || discordantReads.contains(mateRead) || mateRead.isUnmapped())
                continue;

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

        // now all possible discordant and junction mate reads have been collected, so test for overlaps with the min/max aligned position
        // process in order of closest to furthest-out reads in the ref base direction
        List<NonJunctionRead> sortedCandidateReads = candidateReads.stream()
                .sorted(Comparator.comparingInt(x -> isForwardJunction ? -x.read().alignmentEnd() : x.read().alignmentStart()))
                .collect(Collectors.toList());

        List<NonJunctionRead> candidateNonJunctionReads = Lists.newArrayList();

        boolean hasGapped = false;

        for(NonJunctionRead njRead : sortedCandidateReads)
        {
            Read read = njRead.read();

            if(isForwardJunction)
            {
                if(read.alignmentEnd() < minAlignedPosition + ASSEMBLY_REF_SIDE_OVERLAP_BASES)
                {
                    hasGapped = true;
                }
                else
                {
                    minAlignedPosition = min(minAlignedPosition, read.alignmentStart());
                }
            }
            else
            {
                if(read.alignmentStart() > maxAlignedPosition - ASSEMBLY_REF_SIDE_OVERLAP_BASES)
                {
                    hasGapped = true;
                }
                else
                {
                    maxAlignedPosition = max(maxAlignedPosition, read.alignmentEnd());
                }
            }

            if(hasGapped)
            {
                if(njRead.type() == CANDIDATE_DISCORDANT)
                    discordantReads.remove(read);
            }
            else
            {
                assembly.checkAddRefSideSoftClip(read);
                candidateNonJunctionReads.add(njRead);
            }
        }

        // for now add these candidate discordant reads without any further checks
        for(NonJunctionRead njRead : candidateNonJunctionReads)
        {
            assembly.addCandidateSupport(njRead.read(), njRead.type());
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

            if(read.alignmentEnd() >= junctionPosition)
                return false;
        }
        else
        {
            if(read.positiveStrand())
                return false;

            if(read.alignmentStart() <= junctionPosition)
                return false;
        }

        return true;
    }

    private boolean isDiscordantCandidate(final Read read, boolean isForwardJunction, int junctionPosition)
    {
        return isValidSupportCoordsVsJunction(read, isForwardJunction, junctionPosition) && isDiscordantFragment(read);
    }

    public static void extendRefBases(
            final JunctionAssembly assembly, final List<SupportRead> nonJunctionSupport, final RefGenomeInterface refGenome,
            boolean allowBranching)
    {
        // find the maximal ref base extension point and make note of any recurrent soft-clip points including possible branched assemblies
        int minAlignedPosition = assembly.minAlignedPosition();
        int maxAlignedPosition = assembly.maxAlignedPosition();

        boolean isForwardJunction = assembly.junction().isForward();
        int initialRefPosition = isForwardJunction ? minAlignedPosition : maxAlignedPosition;

        Collections.sort(nonJunctionSupport,
                Comparator.comparingInt(x -> isForwardJunction ? -x.cachedRead().alignmentEnd() : x.cachedRead().alignmentStart()));

        boolean hasGapped = false; // currently unused since gaps are filled in with ref bases

        List<RefSideSoftClip> refSideSoftClips = assembly.refSideSoftClips();
        refSideSoftClips.clear(); // will be re-established

        for(SupportRead support : nonJunctionSupport)
        {
            if(isForwardJunction)
            {
                hasGapped |= support.alignmentEnd() < minAlignedPosition;
                minAlignedPosition = min(minAlignedPosition, support.alignmentStart());
            }
            else
            {
                hasGapped |= support.alignmentStart() > maxAlignedPosition;
                maxAlignedPosition = max(maxAlignedPosition, support.alignmentEnd());
            }

            RefSideSoftClip.checkAddRefSideSoftClip(refSideSoftClips, assembly.junction(), support.cachedRead());
        }

        int nonSoftClipRefPosition = isForwardJunction ? minAlignedPosition : maxAlignedPosition;
        purgeRefSideSoftClips(refSideSoftClips, PRIMARY_ASSEMBLY_MIN_READ_SUPPORT, REF_SIDE_MIN_SOFT_CLIP_LENGTH, nonSoftClipRefPosition);

        boolean hasUnmatched = assembly.refSideSoftClips().stream().anyMatch(x -> !x.matchesOriginal());

        if(!hasUnmatched || !allowBranching)
        {
            if(nonSoftClipRefPosition != initialRefPosition && !nonJunctionSupport.isEmpty())
            {
                RefBaseAssembly refBaseAssembly = new RefBaseAssembly(assembly, nonSoftClipRefPosition, refGenome);

                checkAddRefAssemblySupport(refBaseAssembly, nonJunctionSupport, Collections.emptySet());

                if(refBaseAssembly.supportCount() > 0)
                    assembly.mergeRefBaseAssembly(refBaseAssembly);
            }

            return;
        }

        // now build out any distinct, branching assemblies from ref-base soft-clips
        List<Set<String>> excludedReadIdsList = allocateExcludedReads(assembly, nonJunctionSupport);
        List<SupportRead> initialSupport = Lists.newArrayList(assembly.support());

        List<JunctionAssembly> branchedAssemblies = Lists.newArrayList(assembly);

        // the original assembly is handled first
        for(int i = 0; i < excludedReadIdsList.size(); ++i)
        {
            Set<String> excludedReads = excludedReadIdsList.get(i);

            JunctionAssembly junctionAssembly = null;

            int extensionRefPosition;

            if(i == 0)
            {
                junctionAssembly = assembly;

                // first remove support from the main assembly
                assembly.removeSupportReads(excludedReads);

                // CHECK: can the main assembly really be reduced and removed? or filter later on
                /*
                if(assembly.supportCount() < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
                {
                    branchedAssemblies.clear();
                    continue;
                }
                */

                extensionRefPosition = isForwardJunction ? minAlignedPosition : maxAlignedPosition;
            }
            else
            {
                RefSideSoftClip refSideSoftClip = refSideSoftClips.get(i - 1);

                if(refSideSoftClip.matchesOriginal() || refSideSoftClip.hasProximateMatch(assembly.refBasePosition()))
                    continue;

                junctionAssembly = new JunctionAssembly(assembly, refSideSoftClip, initialSupport, excludedReads);

                if(junctionAssembly.supportCount() < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
                    continue;

                extensionRefPosition = refSideSoftClip.Position;
            }

            RefBaseAssembly refBaseAssembly = new RefBaseAssembly(junctionAssembly, extensionRefPosition, refGenome);

            checkAddRefAssemblySupport(refBaseAssembly, nonJunctionSupport, excludedReads);

            // only add branched assemblies if they have sufficient support
            if(junctionAssembly != assembly)
            {
                if(refBaseAssembly.supportCount() == 0)
                    continue;

                junctionAssembly.mergeRefBaseAssembly(refBaseAssembly);

                if(AssemblyUtils.hasUnsetBases(junctionAssembly))
                    continue;

                // same criteria as above
                if(junctionAssembly.supportCount() < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
                    continue;

                branchedAssemblies.add(junctionAssembly);
                junctionAssembly.setOutcome(DUP_BRANCHED);
            }
            else
            {
                if(refBaseAssembly.supportCount() > 0) // same criteria as above
                {
                    junctionAssembly.mergeRefBaseAssembly(refBaseAssembly);
                }
            }
        }

        // set references between them - for now just for TSV output
        for(JunctionAssembly junctionAssembly : branchedAssemblies)
        {
            if(junctionAssembly != assembly)
                assembly.phaseGroup().addDerivedAssembly(junctionAssembly);
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

    private static void checkAddRefAssemblySupport(
            final RefBaseAssembly refBaseAssembly, final List<SupportRead> nonJunctionReads, final Set<String> excludedReadIds)
    {
        for(SupportRead support : nonJunctionReads)
        {
            if(!excludedReadIds.contains(support.id()))
            {
                SupportType type = support.type() == CANDIDATE_DISCORDANT ? DISCORDANT : support.type();
                refBaseAssembly.checkAddRead(support.cachedRead(), type, ASSEMBLY_EXTENSION_BASE_MISMATCH, ASSEMBLY_REF_SIDE_OVERLAP_BASES);
            }
        }
    }

}
