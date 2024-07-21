package com.hartwig.hmftools.esvee.assembly;

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
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.basesMatch;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.INVALID_INDEX;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.DUP_BRANCHED;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.findIndelExtensions;
import static com.hartwig.hmftools.esvee.assembly.RemoteRegionFinder.findRemoteRegions;
import static com.hartwig.hmftools.esvee.assembly.types.ReadAssemblyIndices.getRefReadIndices;
import static com.hartwig.hmftools.esvee.assembly.types.RefSideSoftClip.purgeRefSideSoftClips;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.DISCORDANT;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION_MATE;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.isDiscordantFragment;
import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.assembly.types.ReadAssemblyIndices;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
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

        int newRefBasePosition = assembly.refBasePosition();

        boolean isForwardJunction = assembly.junction().isForward();
        int junctionPosition = assembly.junction().Position;

        // not keeping reads with unmapped mates since not sure how to incorporate their bases
        List<Read> discordantReads = unfilteredNonJunctionReads.stream()
                .filter(x -> isDiscordantCandidate(x, isForwardJunction, junctionPosition) || x.isMateUnmapped())
                .filter(x -> !assembly.hasReadSupport(x.mateRead()))
                .collect(Collectors.toList());

        List<NonJunctionRead> candidateReads = discordantReads.stream()
                .map(x -> new NonJunctionRead(x, DISCORDANT)).collect(Collectors.toList());

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

            mateRead.markJunctionMate();

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
            assembly.addCandidateSupport(read.read());
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
                if(read.alignmentEnd() < newRefBasePosition + ASSEMBLY_REF_SIDE_OVERLAP_BASES)
                    break;

                newRefBasePosition = min(newRefBasePosition, read.alignmentStart());
            }
            else
            {
                if(read.alignmentStart() > newRefBasePosition - ASSEMBLY_REF_SIDE_OVERLAP_BASES)
                    break;

                newRefBasePosition = max(newRefBasePosition, read.alignmentEnd());
            }

            assembly.checkAddRefSideSoftClip(read);
        }

        findRemoteRegions(assembly, discordantReads, remoteJunctionMates, suppJunctionReads);

        // only keep possible alternative ref-base assemblies with sufficient evidence and length
        purgeRefSideSoftClips(assembly.refSideSoftClips(), PRIMARY_ASSEMBLY_MIN_READ_SUPPORT, REF_SIDE_MIN_SOFT_CLIP_LENGTH, newRefBasePosition);
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
            final JunctionAssembly assembly, final List<Read> candidateSupport, final RefGenomeInterface refGenome, boolean allowBranching)
    {
        if(candidateSupport.isEmpty())
            return;

        // find the maximal ref base extension point and make note of any recurrent soft-clip points including possible branched assemblies

        boolean isForwardJunction = assembly.junction().isForward();
        int initialRefPosition = assembly.refBasePosition();
        int newRefBasePosition = initialRefPosition;

        // build out the reads starting with those closest to the junction, and test for gaps in ref bases
        Collections.sort(candidateSupport,
                Comparator.comparingInt(x -> isForwardJunction ? -x.alignmentEnd() : x.alignmentStart()));

        // capture RSSC from these new candidate reads
        List<RefSideSoftClip> refSideSoftClips = assembly.refSideSoftClips();
        List<Read> nonJunctionSupport = Lists.newArrayListWithExpectedSize(candidateSupport.size());

        for(Read read : candidateSupport)
        {
            if(isForwardJunction)
            {
                newRefBasePosition = min(newRefBasePosition, read.alignmentStart());
            }
            else
            {
                newRefBasePosition = max(newRefBasePosition, read.alignmentEnd());
            }

            nonJunctionSupport.add(read);

            RefSideSoftClip.checkAddRefSideSoftClip(refSideSoftClips, assembly.junction(), read);
        }

        if(nonJunctionSupport.isEmpty())
            return;

        int nonSoftClipRefPosition = newRefBasePosition;

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
            final JunctionAssembly assembly, int newRefPosition, final List<Read> nonJunctionSupport,
            final RefGenomeInterface refGenome, boolean isSoftClipped)
    {
        if(newRefPosition != assembly.refBasePosition())
        {
            if(assembly.isForwardJunction() == (newRefPosition < assembly.refBasePosition()))
            {
                assembly.extendRefBases(newRefPosition, Collections.emptyList(), refGenome);
            }
            else if(isSoftClipped)
            {
                // need to trim the existing ref bases due to ref-side soft-clipped determining the bounds
                assembly.trimRefBases(newRefPosition);
            }
        }

        // test and add support
        checkAddRefBaseSupport(assembly, nonJunctionSupport, Collections.emptySet());
    }

    private static void branchAssembliesFromRefBases(
            final JunctionAssembly originalAssembly, int nonSoftClipRefPosition,
            final List<Read> nonJunctionSupport, int maxRefSideSupport, final RefGenomeInterface refGenome)
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

            if(i == 0)
            {
                junctionAssembly = originalAssembly;

                // first remove support from the main assembly
                originalAssembly.removeSupportReads(excludedReads);

                extendAssemblyRefBases(originalAssembly, nonSoftClipRefPosition, nonJunctionSupport, refGenome, false);
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

                List<SupportRead> newSupport = initialSupport.stream().filter(x -> !excludedReads.contains(x.id())).collect(Collectors.toList());

                junctionAssembly = new JunctionAssembly(originalAssembly, refSideSoftClip, newRefBaseLength, newSupport);
                extendAssemblyRefBases(junctionAssembly, refSideSoftClip.Position, nonJunctionSupport, refGenome, true);
            }

            checkAddRefBaseSupport(junctionAssembly, nonJunctionSupport, excludedReads);

            // only add branched assemblies if they have sufficient support
            if(junctionAssembly != originalAssembly)
            {
                // check if has sufficient support to branch the assembly
                int totalSupport = junctionAssembly.supportCount();
                double supportPerc = totalSupport / (double)maxRefSideSupport;
                boolean hasSufficientSecondRefSupport = totalSupport >= PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT
                        && supportPerc >= PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC;

                if(!hasSufficientSecondRefSupport)
                    continue;

                // junctionAssembly.buildRepeatInfo();
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

    private static void checkAddRefBaseSupport(
            final JunctionAssembly assembly, final List<Read> nonJunctionReads, final Set<String> excludedReadIds)
    {
        List<Read> newSupportReads = !excludedReadIds.isEmpty() ?
                nonJunctionReads.stream().filter(x -> !excludedReadIds.contains(x.id())).collect(Collectors.toList()) : nonJunctionReads;

        // favour junction mates first, then reads with least variants and most aligned bases
        Collections.sort(newSupportReads, new RefBaseReadComparator());

        for(Read read : newSupportReads)
        {
            SupportType type = read.hasJunctionMate() ? JUNCTION_MATE : DISCORDANT;
            checkAddRefBaseRead(assembly, read, type);
        }
    }

    private static final int REF_READ_SEARCH_LENGTH = 20;

    public static boolean checkAddRefBaseRead(final JunctionAssembly assembly, final Read read, final SupportType supportType)
    {
        ReadAssemblyIndices readIndexInfo = getRefReadIndices(assembly, assembly.refBasePosition(), read);

        if(readIndexInfo == null)
            return false;

        int readStartIndex = readIndexInfo.ReadIndexStart;
        final byte[] assemblyBases = assembly.bases();
        final byte[] assemblyBaseQuals = assembly.baseQuals();

        boolean canAddRead = canAddRefBaseRead(assemblyBases, assemblyBaseQuals, read, readIndexInfo);

        if(!canAddRead && readStartIndex < read.getBases().length - REF_READ_SEARCH_LENGTH)
        {
            // run a simple sequence search to find the alignment start where prior indels have offset the read's infer assembly index start
            int readTestEndIndex = readStartIndex + REF_READ_SEARCH_LENGTH;
            int length = readTestEndIndex - readStartIndex + 1;

            if(readStartIndex < 0 || readStartIndex >= read.getBases().length || readStartIndex + length > read.getBases().length)
            {
                SV_LOGGER.error("refAssembly({}) invalid indices({} - {}) vs readBases({}) for ref extension read search",
                        assembly, readStartIndex, readTestEndIndex, read.getBases().length);

                System.exit(1);
            }

            String readBases = new String(read.getBases(), readStartIndex, length);
            int assemblyStartIndex = new String(assembly.bases()).indexOf(readBases);

            if(assemblyStartIndex >= 0)
            {
                readIndexInfo = new ReadAssemblyIndices(readStartIndex, readIndexInfo.ReadIndexEnd, assemblyStartIndex);
                canAddRead = canAddRefBaseRead(assemblyBases, assemblyBaseQuals, read, readIndexInfo);
            }
        }

        if(!canAddRead)
        {
            // junction mate reads are added as support even if their ref bases don't match
            if(supportType == SupportType.JUNCTION_MATE)
            {
                SupportRead supportRead = new SupportRead(
                        read, supportType, INVALID_INDEX, 0, ASSEMBLY_EXTENSION_BASE_MISMATCH + 1);
                assembly.support().add(supportRead);
            }

            return false;
        }

        assembly.addRead(read, readIndexInfo, supportType, null);

        return true;
    }

    private static boolean canAddRefBaseRead(
            final byte[] assemblyBases, final byte[] assemblyBaseQuals, final Read read, final ReadAssemblyIndices readIndexInfo)
    {
        int mismatchCount = 0;
        int overlappedBaseCount = 0;
        int assemblyIndex = readIndexInfo.AssemblyIndexStart;

        int permittedMismatches = ASSEMBLY_EXTENSION_BASE_MISMATCH;
        int requiredOverlap = ASSEMBLY_REF_SIDE_OVERLAP_BASES;

        for(int i = readIndexInfo.ReadIndexStart; i <= readIndexInfo.ReadIndexEnd; ++i, ++assemblyIndex)
        {
            if(assemblyIndex >= assemblyBases.length)
                break;

            if(assemblyBases[assemblyIndex] == 0)
                continue;

            ++overlappedBaseCount;

            // any unset base (ie unset qual) can be a mismatch
            byte refBaseQual = assemblyBaseQuals[assemblyIndex] == 0 ? (byte)(LOW_BASE_QUAL_THRESHOLD + 1) : assemblyBaseQuals[assemblyIndex];

            if(!basesMatch(
                    read.getBases()[i], assemblyBases[assemblyIndex], read.getBaseQuality()[i], refBaseQual, LOW_BASE_QUAL_THRESHOLD))
            {
                ++mismatchCount;

                if(mismatchCount > permittedMismatches)
                    return false;
            }
        }

        return overlappedBaseCount >= requiredOverlap;
    }

    private static class RefBaseReadComparator implements Comparator<Read>
    {
        @Override
        public int compare(Read first, Read second)
        {
            if(first.hasJunctionMate() != second.hasJunctionMate())
                return first.hasJunctionMate() ? -1 : 1;

            // favour reads with less variants
            if(first.numOfEvents() != second.numOfEvents())
                return first.numOfEvents() < second.numOfEvents() ? -1 : 1;

            int firstAlignedLength = first.alignmentEnd() - first.alignmentStart() + 1;
            int secondAlignedLength = second.alignmentEnd() - second.alignmentStart() + 1;
            return -1 * Integer.compare(firstAlignedLength, secondAlignedLength);
        }
    }

    private static List<Set<String>> allocateExcludedReads(final JunctionAssembly assembly, final List<Read> nonJunctionReads)
    {
        // split up reads into those supporting the new branched assemblies vs the original
        List<RefSideSoftClip> refSideSoftClips = assembly.refSideSoftClips();

        boolean isForwardJunction = assembly.isForwardJunction();
        Set<String> nscReads = Sets.newHashSet();

        int maxRefSideSoftClipPosition = isForwardJunction ?
                refSideSoftClips.stream().filter(x -> !x.matchesOriginal()).mapToInt(x -> x.Position).min().orElse(0)
                : refSideSoftClips.stream().filter(x -> !x.matchesOriginal()).mapToInt(x -> x.Position).max().orElse(0);

        Set<String> softClippedReadIds = Sets.newHashSet();
        refSideSoftClips.stream().filter(x -> !x.matchesOriginal()).forEach(x -> softClippedReadIds.addAll(x.readIds()));

        // find any supporting read or read mate which extends past the furthest soft-clip position
        for(SupportRead read : assembly.support())
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

        for(Read read : nonJunctionReads)
        {
            if(softClippedReadIds.stream().anyMatch(x -> x.equals(read.id())))
                continue;

            boolean addRead = false;

            if(isForwardJunction)
            {
                if(read.alignmentStart() < maxRefSideSoftClipPosition
                || (read.hasJunctionMate() && read.isMateMapped() && read.mateAlignmentStart() < maxRefSideSoftClipPosition))
                {
                    addRead = true;
                }
            }
            else
            {
                if(read.alignmentEnd() > maxRefSideSoftClipPosition
                || (read.hasJunctionMate() && read.isMateMapped() && read.mateAlignmentEnd() > maxRefSideSoftClipPosition))
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

    public static void trimAssemblyRefBases(final JunctionAssembly assembly, final int maxGap)
    {
        // trim back any ref bases where the gap between reads exceeds the permitted length, and any unset bases
        final byte[] assemblyBases = assembly.bases();
        final byte[] assemblyBaseQuals = assembly.baseQuals();

        int refBasePosition = assembly.refBasePosition();
        int newRefBasePosition = refBasePosition;

        if(assembly.isForwardJunction())
        {
            int trimIndex = 0;
            int lastSetBase = -1;

            for(int i = 0; i < assembly.junctionIndex(); ++i)
            {
                if(assemblyBaseQuals[i] > 0)
                {
                    if(lastSetBase >= 0 && i - lastSetBase > maxGap)
                        trimIndex = i;

                    lastSetBase = i;
                }
                else if(lastSetBase < 0)
                {
                    trimIndex = i + 1;
                }
            }

            newRefBasePosition += trimIndex;
        }
        else
        {
            int lastIndex = assemblyBases.length - 1;;
            int trimIndex = lastIndex;
            int lastSetBase = -1;

            for(int i = lastIndex; i > assembly.junctionIndex(); --i)
            {
                if(assemblyBaseQuals[i] > 0)
                {
                    if(lastSetBase >= 0 && lastSetBase - i > maxGap)
                        trimIndex = i;

                    lastSetBase = i;
                }
                else if(lastSetBase < 0)
                {
                    trimIndex = i - 1;
                }
            }

            newRefBasePosition -= (lastIndex - trimIndex);
        }

        if(newRefBasePosition != refBasePosition)
            assembly.trimRefBases(newRefBasePosition);
    }
}
