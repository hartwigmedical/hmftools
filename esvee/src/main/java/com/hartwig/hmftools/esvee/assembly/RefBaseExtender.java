package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.inferredInsertSizeAbs;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_EXTENSION_BASE_MISMATCH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_REF_BASE_MAX_GAP;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_READ_OVERLAP_BASES;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_SPLIT_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.UNMAPPED_TRIM_THRESHOLD;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.basesMatch;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.isValidSupportCoordsVsJunction;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.DUP_BRANCHED;
import static com.hartwig.hmftools.esvee.assembly.RemoteRegionFinder.findRemoteRegions;
import static com.hartwig.hmftools.esvee.assembly.types.ReadAssemblyIndices.getRefReadIndices;
import static com.hartwig.hmftools.esvee.assembly.types.RefSideSoftClip.checkAddRefSideSoftClip;
import static com.hartwig.hmftools.esvee.assembly.types.RefSideSoftClip.purgeRefSideSoftClips;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.DISCORDANT;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION_MATE;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.isDiscordantFragment;
import static com.hartwig.hmftools.esvee.common.CommonUtils.aboveMinQual;
import static com.hartwig.hmftools.esvee.common.CommonUtils.isDuplicationFragment;
import static com.hartwig.hmftools.esvee.common.SvConstants.DEFAULT_MAX_CONCORDANT_FRAG_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.assembly.read.ReadAdjustments;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.ReadAssemblyIndices;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.RefSideSoftClip;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;
import com.hartwig.hmftools.esvee.assembly.read.Read;

public class RefBaseExtender
{
    public RefBaseExtender() { }

    public static void checkRefSideSoftClips(final JunctionAssembly assembly)
    {
        // purge any junction support which extends beyond a consensus ref-side soft-clip
        List<RefSideSoftClip> refSideSoftClips = Lists.newArrayList();

        for(SupportRead read : assembly.support())
        {
            checkAddRefSideSoftClip(refSideSoftClips, assembly.junction(), read.cachedRead());
        }

        if(refSideSoftClips.isEmpty())
            return;

        Collections.sort(refSideSoftClips, Comparator.comparingInt(x -> -x.readCount()));

        RefSideSoftClip refSideSoftClip = refSideSoftClips.get(0);

        if(refSideSoftClip.readCount() < ASSEMBLY_MIN_READ_SUPPORT)
            return;

        int mainSoftClipCount = refSideSoftClip.readCount();
        int totalSoftClipCount = refSideSoftClips.stream().mapToInt(x -> x.readCount()).sum();
        int nonSoftClipCount = assembly.supportCount() - totalSoftClipCount;

        if(nonSoftClipCount > 0 && nonSoftClipCount < mainSoftClipCount && nonSoftClipCount < ASSEMBLY_SPLIT_MIN_READ_SUPPORT)
        {
            // as per the branching routine run during linking, require a minimum number of reads to keep both reads which soft-clip and
            // those which run past that point
            List<SupportRead> support = assembly.support();
            int index = 0;

            while(index < support.size())
            {
                SupportRead read = support.get(index);

                if(refSideSoftClips.stream().noneMatch(x -> x.readIds().contains(read.id())))
                {
                    support.remove(index);
                }
                else
                {
                    ++index;
                }
            }
        }
    }

    public void findAssemblyCandidateExtensions(final JunctionAssembly assembly, final List<Read> unfilteredNonJunctionReads)
    {
        // find all possible discordant reads and junction mate reads, and use them to extend the ref bases
        // other applicable info such as soft-clips on the ref side and links to remote regions are also captured
        int newRefBasePosition = assembly.refBasePosition();
        boolean isForwardJunction = assembly.junction().isForward();
        int junctionPosition = assembly.junction().Position;

        // difference for local cigar-based indels
        boolean isIndelJunction = assembly.indel();

        List<NonJunctionRead> candidateReads = Lists.newArrayList();
        List<Read> discordantReads;

        if(isIndelJunction)
        {
            discordantReads = Collections.emptyList();
        }
        else
        {
            Set<String> concordantReadIds = unfilteredNonJunctionReads.stream().filter(x -> isConcordantRead(x))
                    .map(x -> x.id()).collect(Collectors.toSet());

            discordantReads = unfilteredNonJunctionReads.stream()
                    .filter(x -> isDiscordantCandidate(x, isForwardJunction, junctionPosition, assembly, concordantReadIds))
                    .collect(Collectors.toList());

            discordantReads.forEach(x -> candidateReads.add(new NonJunctionRead(x, DISCORDANT)));

            discordantReads.stream().filter(x -> x.isMateUnmapped() && x.mateRead() != null && !filterUnmapped(x.mateRead(), true))
                    .forEach(x -> assembly.addUnmappedRead(x.mateRead()));
        }

        List<Read> remoteJunctionMates = Lists.newArrayList();
        List<Read> suppJunctionReads = Lists.newArrayList();

        // add any junction mates in the same window
        for(SupportRead read : assembly.support())
        {
            if(read.cachedRead().hasSupplementary())
                suppJunctionReads.add(read.cachedRead());

            if(!isIndelJunction && read.isDiscordant())
            {
                remoteJunctionMates.add(read.cachedRead());
                continue;
            }

            // look to extend from local mates on the ref side of the junction
            Read mateRead = read.cachedRead().mateRead();

            if(mateRead == null)
                continue;

            mateRead.markJunctionMate();

            if(discordantReads.contains(mateRead))
                continue;

            if(!mateRead.isUnmapped())
            {
                boolean isPastJunction = (isForwardJunction && mateRead.alignmentEnd() >= junctionPosition)
                        || (!isForwardJunction && mateRead.alignmentStart() <= junctionPosition);

                if(isPastJunction)
                {
                    if(!mateRead.isUnmapped() && !mateRead.isLeftClipped() && !mateRead.isRightClipped())
                        assembly.addConcordantCandidate(mateRead);
                }
                else
                {
                    candidateReads.add(new NonJunctionRead(mateRead, JUNCTION_MATE));
                }
            }
            else
            {
                if(filterUnmapped(mateRead,false))
                    continue;

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

            boolean exceedsDiscGap = false;

            if(isForwardJunction)
            {
                if(!read.hasJunctionMate() && read.alignmentEnd() < newRefBasePosition - ASSEMBLY_REF_BASE_MAX_GAP)
                    exceedsDiscGap = true;
                else
                    newRefBasePosition = min(newRefBasePosition, read.alignmentStart());
            }
            else
            {
                if(!read.hasJunctionMate() && read.alignmentStart() > newRefBasePosition + ASSEMBLY_REF_BASE_MAX_GAP)
                    exceedsDiscGap = true;
                else
                    newRefBasePosition = max(newRefBasePosition, read.alignmentEnd());
            }

            if(!exceedsDiscGap)
                assembly.checkAddRefSideSoftClip(read);
        }

        // consolidate all links to remote regions for later use in phase group building and assembly linking
        if(!isIndelJunction)
            findRemoteRegions(assembly, discordantReads, remoteJunctionMates, suppJunctionReads);

        // only keep possible alternative ref-base assemblies with sufficient evidence and length
        purgeRefSideSoftClips(assembly.refSideSoftClips(), newRefBasePosition);
    }

    private static boolean isConcordantRead(final Read read)
    {
        if(!read.isPairedRead() || !read.isMateMapped() || !read.isMateMapped())
            return false;

        if(!read.chromosome().equals(read.mateChromosome()) || read.orientation() == read.mateOrientation())
            return false;

        int fragmentSize = inferredInsertSizeAbs(read.bamRecord());

        if(fragmentSize > DEFAULT_MAX_CONCORDANT_FRAG_LENGTH)
            return false;

        if(isDuplicationFragment(read.bamRecord(), fragmentSize))
            return false;

        return true;
    }

    private static boolean isDiscordantCandidate(
            final Read read, boolean isForwardJunction, int junctionPosition, final JunctionAssembly assembly, final Set<String> concordantReadIds)
    {
        if(!isValidSupportCoordsVsJunction(read, isForwardJunction, junctionPosition))
            return false;

        if(!isDiscordantFragment(read))
        {
            if(!read.isMateUnmapped())
                return false;

            // test the mate read's base quals and
            if(read.mateRead() != null && filterUnmapped(read.mateRead(), true))
                return false;
        }
        else
        {
            // must not match a concordant fragment as determined by a local supplementary or vice-versa
            if(concordantReadIds.contains(read.id()))
                return false;
        }

        // skip if the read's mate is a junction read
        return !assembly.hasReadSupport(read.mateRead());
    }

    private static boolean filterUnmapped(final Read read, boolean isDiscordant)
    {
        if(!read.isUnmapped())
            return false;

        read.trimLowQualBases();

        if(isDiscordant && read.mateRead() != null && read.mateRead().mappingQuality() == 0)
            return true;

        return read.basesLength() - read.baseTrimCount() < UNMAPPED_TRIM_THRESHOLD;
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

    public static void extendRefBases(
            final JunctionAssembly assembly, final List<Read> candidateSupport, final RefGenomeInterface refGenome, boolean allowBranching,
            boolean allowSoftClipRestrictions)
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
        // NOTE: this is only done once per assembly linking and extension for now to avoid repeated consideration of branching
        boolean considerRefSideSoftClips = allowSoftClipRestrictions && candidateSupport.stream().anyMatch(x -> x.hasJunctionMate());

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

            if(considerRefSideSoftClips)
                checkAddRefSideSoftClip(refSideSoftClips, assembly.junction(), read);
        }

        if(nonJunctionSupport.isEmpty())
            return;

        int nonSoftClipRefPosition = newRefBasePosition;

        if(considerRefSideSoftClips)
            purgeRefSideSoftClips(refSideSoftClips, nonSoftClipRefPosition);

        if(!considerRefSideSoftClips || refSideSoftClips.isEmpty())
        {
            // most common scenario
            extendAssemblyRefBases(assembly, nonSoftClipRefPosition, nonJunctionSupport, refGenome, false);
            return;
        }

        Collections.sort(refSideSoftClips, Comparator.comparingInt(x -> -x.readCount()));

        RefSideSoftClip candidateRefSideSoftClip = refSideSoftClips.get(0);

        int nonSoftClipSupport = 0;
        int junctionSoftClipped = 0;
        int junctionNonSoftClipped = 0;

        for(SupportRead read : assembly.support())
        {
            if(isForwardJunction)
            {
                if(read.isLeftClipped() && read.alignmentStart() == candidateRefSideSoftClip.Position)
                    ++junctionSoftClipped;
                else if(!read.isLeftClipped() && read.alignmentStart() < candidateRefSideSoftClip.Position)
                    ++junctionNonSoftClipped;
            }
            else
            {
                if(read.isRightClipped() && read.alignmentEnd() == candidateRefSideSoftClip.Position)
                    ++junctionSoftClipped;
                else if(!read.isRightClipped() && read.alignmentEnd() > candidateRefSideSoftClip.Position)
                    ++junctionNonSoftClipped;
            }
        }

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

        boolean majorityJuncReadsSoftClipped = junctionSoftClipped > junctionNonSoftClipped;

        if(majorityJuncReadsSoftClipped || nonSoftClipSupport <= candidateRefSideSoftClip.readCount())
        {
            primaryRefPosition = candidateRefSideSoftClip.Position;
            primaryRefPositionSupport = candidateRefSideSoftClip.readCount();
            secondRefPositionSupport = nonSoftClipSupport;
            usesSoftClippedPosition = true;
        }
        else
        {
            primaryRefPosition = nonSoftClipRefPosition;
            primaryRefPositionSupport = nonSoftClipSupport;
            secondRefPositionSupport = candidateRefSideSoftClip.readCount();
        }

        double secondRefPositionSupportPerc = secondRefPositionSupport / (double)primaryRefPositionSupport;
        boolean hasSufficientSecondRefSupport = secondRefPositionSupport >= ASSEMBLY_SPLIT_MIN_READ_SUPPORT
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
                assembly.extendRefBases(newRefPosition, refGenome);
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

            List<Read> nonExcludedNonJunctionSupport = nonJunctionSupport.stream()
                    .filter(x -> excludedReads.stream().noneMatch(y -> x.id().contains(y))).collect(Collectors.toList());

            JunctionAssembly junctionAssembly = null;

            if(i == 0)
            {
                junctionAssembly = originalAssembly;

                // first remove support from the main assembly
                originalAssembly.removeSupportReads(excludedReads);

                extendAssemblyRefBases(originalAssembly, nonSoftClipRefPosition, nonExcludedNonJunctionSupport, refGenome, false);
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

                extendAssemblyRefBases(junctionAssembly, refSideSoftClip.Position, nonExcludedNonJunctionSupport, refGenome, true);
            }

            // only add branched assemblies if they have sufficient support
            if(junctionAssembly != originalAssembly)
            {
                // check if has sufficient support to branch the assembly
                int totalSupport = junctionAssembly.supportCount();
                double supportPerc = totalSupport / (double)maxRefSideSupport;
                boolean hasSufficientSecondRefSupport = totalSupport >= ASSEMBLY_SPLIT_MIN_READ_SUPPORT
                        && supportPerc >= PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC;

                if(!hasSufficientSecondRefSupport)
                    continue;

                branchedAssemblies.add(junctionAssembly);
                junctionAssembly.setOutcome(DUP_BRANCHED);

                // for now only 1 branched assembly will be made
                break;
            }
        }

        // set references between them - for now just for TSV output
        for(JunctionAssembly junctionAssembly : branchedAssemblies)
        {
            purgeUnrelatedRefSideSoftClips(junctionAssembly);

            if(junctionAssembly != originalAssembly)
                originalAssembly.phaseGroup().addDerivedAssembly(junctionAssembly);
        }
    }

    private static void purgeUnrelatedRefSideSoftClips(final JunctionAssembly assembly)
    {
        int index = 0;
        while(index < assembly.refSideSoftClips().size())
        {
            RefSideSoftClip refSideSoftClip = assembly.refSideSoftClips().get(index);

            if(refSideSoftClip.hasProximateMatch(assembly.refBasePosition()))
                ++index;
            else
                assembly.refSideSoftClips().remove(index);
        }
    }

    private static void checkAddRefBaseSupport(
            final JunctionAssembly assembly, final List<Read> nonJunctionReads, final Set<String> excludedReadIds)
    {
        List<Read> secondarySupportReads = Lists.newArrayList();

        int permittedMismatches = ASSEMBLY_EXTENSION_BASE_MISMATCH;

        for(Read read : nonJunctionReads)
        {
            if(excludedReadIds.contains(read.id()))
                continue;

            // favour junction mates with ref base support first
            if(read.hasJunctionMate())
            {
                SupportRead juncMate = assembly.support().stream()
                        .filter(x -> x.type() == JUNCTION && x.matchesFragment(read, false)).findFirst().orElse(null);

                if(juncMate != null && juncMate.hasReferenceMismatches() && juncMate.referenceMismatches() <= permittedMismatches)
                {
                    checkAddRefBaseRead(assembly, read, JUNCTION_MATE, 0);
                    continue;
                }
            }

            secondarySupportReads.add(read);
        }

        boolean sortOnReadStart = assembly.isReverseJunction();
        Collections.sort(secondarySupportReads, Comparator.comparingInt(x -> sortOnReadStart ? x.alignmentStart() : -x.alignmentEnd()));

        for(Read read : secondarySupportReads)
        {
            SupportType type = read.hasJunctionMate() ? JUNCTION_MATE : DISCORDANT;
            checkAddRefBaseRead(assembly, read, type);
        }
    }

    private static final int REF_READ_SEARCH_LENGTH = 20;

    public static boolean checkAddRefBaseRead(final JunctionAssembly assembly, final Read read, final SupportType supportType)
    {
        return checkAddRefBaseRead(assembly, read, supportType, ASSEMBLY_READ_OVERLAP_BASES);
    }

    private static boolean checkAddRefBaseRead(
            final JunctionAssembly assembly, final Read read, final SupportType supportType, int requiredOverlap)
    {
        ReadAssemblyIndices readAssemblyIndices = getRefReadIndices(assembly, assembly.refBasePosition(), read);

        if(readAssemblyIndices == null)
            return false;

        int readStartIndex = readAssemblyIndices.ReadIndexStart;
        final byte[] assemblyBases = assembly.bases();
        final byte[] assemblyBaseQuals = assembly.baseQuals();

        int permittedMismatches = ASSEMBLY_EXTENSION_BASE_MISMATCH;

        boolean canAddRead = canAddRefBaseRead(assemblyBases, assemblyBaseQuals, read, readAssemblyIndices, requiredOverlap, permittedMismatches);

        if(!canAddRead && readStartIndex < read.getBases().length - REF_READ_SEARCH_LENGTH)
        {
            // run a simple sequence search to find the alignment start where prior indels have offset the read's infer assembly index start
            int readTestEndIndex = readStartIndex + REF_READ_SEARCH_LENGTH;
            int length = readTestEndIndex - readStartIndex + 1;

            if(readStartIndex < 0 || readStartIndex >= read.getBases().length || readStartIndex + length > read.getBases().length)
            {
                SV_LOGGER.error("refAssembly({}) invalid indices({} - {}) vs readBases({}) for ref extension read search",
                        assembly, readStartIndex, readTestEndIndex, read.getBases().length);
                return false;
            }

            String readBases = new String(read.getBases(), readStartIndex, length);
            int assemblyStartIndex = new String(assembly.bases()).indexOf(readBases);

            if(assemblyStartIndex >= 0)
            {
                readAssemblyIndices = new ReadAssemblyIndices(readStartIndex, readAssemblyIndices.ReadIndexEnd, assemblyStartIndex);

                canAddRead = canAddRefBaseRead(
                        assemblyBases, assemblyBaseQuals, read, readAssemblyIndices, requiredOverlap, permittedMismatches);
            }
        }

        if(!canAddRead)
        {
            // junction mate reads are added as support even if their ref bases don't match
            if(supportType == SupportType.JUNCTION_MATE)
            {
                int junctionReadStartDistance = readAssemblyIndices.junctionReadStartDistance(assembly.junctionIndex());

                SupportRead supportRead = new SupportRead(
                        read, supportType, junctionReadStartDistance, 0, ASSEMBLY_EXTENSION_BASE_MISMATCH + 1);
                assembly.support().add(supportRead);
            }

            return false;
        }

        assembly.addRead(read, readAssemblyIndices, supportType);

        return true;
    }

    private static boolean canAddRefBaseRead(
            final byte[] assemblyBases, final byte[] assemblyBaseQuals, final Read read, final ReadAssemblyIndices readIndexInfo,
            int requiredOverlap, int permittedMismatches)
    {
        int mismatchCount = 0;
        int overlappedBaseCount = 0;
        int assemblyIndex = readIndexInfo.AssemblyIndexStart;

        for(int i = readIndexInfo.ReadIndexStart; i <= readIndexInfo.ReadIndexEnd; ++i, ++assemblyIndex)
        {
            if(i < 0 || assemblyIndex >= assemblyBases.length || i >= read.getBases().length)
                break;

            if(assemblyBases[assemblyIndex] == 0)
                continue;

            if(aboveMinQual(assemblyBaseQuals[assemblyIndex]))
                ++overlappedBaseCount;

            // any unset base (ie unset qual) can be a mismatch
            byte refBaseQual = assemblyBaseQuals[assemblyIndex] == 0 ? (byte)(LOW_BASE_QUAL_THRESHOLD + 1) : assemblyBaseQuals[assemblyIndex];

            if(!basesMatch(read.getBases()[i], assemblyBases[assemblyIndex], read.getBaseQuality()[i], refBaseQual))
            {
                ++mismatchCount;

                if(mismatchCount > permittedMismatches)
                    return false;
            }
        }

        return overlappedBaseCount >= requiredOverlap;
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
