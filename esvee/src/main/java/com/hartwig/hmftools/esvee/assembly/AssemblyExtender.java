package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.getMateAlignmentEnd;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.ASSEMBLY_EXTENSION_BASE_MISMATCH;
import static com.hartwig.hmftools.esvee.SvConstants.ASSEMBLY_EXTENSION_OVERLAP_BASES;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_LENGTH;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.common.AssemblyUtils.findUnsetBases;
import static com.hartwig.hmftools.esvee.common.RefSideSoftClip.purgeRefSideSoftClips;
import static com.hartwig.hmftools.esvee.common.RemoteRegion.REMOTE_READ_TYPE_DISCORDANT_READ;
import static com.hartwig.hmftools.esvee.common.RemoteRegion.REMOTE_READ_TYPE_JUNCTION_MATE;
import static com.hartwig.hmftools.esvee.common.RemoteRegion.REMOTE_READ_TYPE_JUNCTION_SUPP;
import static com.hartwig.hmftools.esvee.common.RemoteRegion.mergeRegions;
import static com.hartwig.hmftools.esvee.common.RemoteRegion.purgeWeakSuppRegions;
import static com.hartwig.hmftools.esvee.common.SupportType.DISCORDANT;
import static com.hartwig.hmftools.esvee.common.SupportType.JUNCTION_MATE;
import static com.hartwig.hmftools.esvee.read.ReadUtils.isDiscordant;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.esvee.common.AssemblySupport;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.RefBaseAssembly;
import com.hartwig.hmftools.esvee.common.RefSideSoftClip;
import com.hartwig.hmftools.esvee.common.RemoteRegion;
import com.hartwig.hmftools.esvee.common.SupportType;
import com.hartwig.hmftools.esvee.read.Read;

public class AssemblyExtender
{
    private final List<JunctionAssembly> mAssemblies;

    public AssemblyExtender(final JunctionAssembly assembly)
    {
        mAssemblies = Lists.newArrayList(assembly);
    }

    public List<JunctionAssembly> assemblies() { return mAssemblies; }

    public void extendAssembly(final List<Read> unfilteredNonJunctionReads)
    {
        // first establish potential boundaries for extending the assembly on the non-junction side
        JunctionAssembly assembly = mAssemblies.get(0);

        int minAlignedPosition = assembly.minAlignedPosition();
        int maxAlignedPosition = assembly.maxAlignedPosition();

        boolean isForwardJunction = assembly.junction().isForward();
        int junctionPosition = assembly.junction().Position;

        List<Read> discordantReads = unfilteredNonJunctionReads.stream()
                .filter(x -> isDiscordantCandidate(x, isForwardJunction, junctionPosition)) // not keeping reads with unmapped mates since not sure how to incorporate their bases
                // .filter(x -> !x.hasMateSet()) // will now be set for local INVs and DUPs, so cannot bea an excluding criteria
                .filter(x -> !assembly.hasReadSupport(x.mateRead()))
                .collect(Collectors.toList());

        List<NonJunctionRead> candidateReads = discordantReads.stream()
                .map(x -> new NonJunctionRead(x, DISCORDANT)).collect(Collectors.toList());

        List<Read> remoteJunctionMates = Lists.newArrayList();
        List<Read> suppJunctionReads = Lists.newArrayList();

        // add any junction mates in the same window
        for(AssemblySupport support : assembly.support())
        {
            if(support.read().hasSupplementary())
                suppJunctionReads.add(support.read());

            if(isDiscordant(support.read()))
            {
                remoteJunctionMates.add(support.read());
                continue;
            }

            // look to extend from local mates on the ref side of the junction
            Read mateRead = support.read().mateRead();

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
                if(read.alignmentEnd() < minAlignedPosition + ASSEMBLY_EXTENSION_OVERLAP_BASES)
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
                if(read.alignmentStart() > maxAlignedPosition - ASSEMBLY_EXTENSION_OVERLAP_BASES)
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
                if(njRead.type() == DISCORDANT)
                    discordantReads.remove(read);
            }
            else
            {
                assembly.checkAddRefSideSoftClip(read);
                candidateNonJunctionReads.add(njRead);
            }
        }

        int nonSoftClipRefPosition = isForwardJunction ? minAlignedPosition : maxAlignedPosition;

        // only keep possible alternative ref-base assemblies with sufficient evidence and length
        List<RefSideSoftClip> refSideSoftClips = assembly.refSideSoftClips();

        purgeRefSideSoftClips(refSideSoftClips, PRIMARY_ASSEMBLY_MIN_READ_SUPPORT, PRIMARY_ASSEMBLY_MIN_LENGTH, nonSoftClipRefPosition);

        boolean hasUnmatched = assembly.refSideSoftClips().stream().anyMatch(x -> !x.matchesOriginal());

        if(!hasUnmatched)
        {
            RefBaseAssembly refBaseAssembly = new RefBaseAssembly(assembly, nonSoftClipRefPosition);

            checkAddRefAssemblySupport(refBaseAssembly, candidateNonJunctionReads, Collections.emptySet());

            if(refBaseAssembly.supportCount() > 0)
                assembly.mergeRefBaseAssembly(refBaseAssembly);

            findRemoteRegions(assembly, Collections.emptySet(), discordantReads, remoteJunctionMates, suppJunctionReads);
            return;
        }

        // now build out any distinct, branching assemblies from ref-base soft-clips
        List<Set<Read>> excludedReadsList = allocateExcludedReads(assembly, candidateNonJunctionReads);
        List<AssemblySupport> initialSupport = Lists.newArrayList(assembly.support());

        // the original assembly is handled first
        for(int i = 0; i < excludedReadsList.size(); ++i)
        {
            Set<Read> excludedReads = excludedReadsList.get(i);

            JunctionAssembly junctionAssembly = null;

            int extensionRefPosition;

            if(i == 0)
            {
                junctionAssembly = assembly;

                // first remove support from the main assembly
                excludedReads.forEach(x -> assembly.removeSupportRead(x));

                // and check if it
                if(assembly.supportCount() < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
                {
                    mAssemblies.clear();
                    continue;
                }

                extensionRefPosition = isForwardJunction ? minAlignedPosition : maxAlignedPosition;
            }
            else
            {
                RefSideSoftClip refSideSoftClip = refSideSoftClips.get(i - 1);

                if(refSideSoftClip.matchesOriginal())
                    continue;

                junctionAssembly = new JunctionAssembly(assembly, refSideSoftClip, initialSupport, excludedReads);

                if(junctionAssembly.supportCount() < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
                    continue;

                extensionRefPosition = refSideSoftClip.Position;
            }

            RefBaseAssembly refBaseAssembly = new RefBaseAssembly(junctionAssembly, extensionRefPosition);

            checkAddRefAssemblySupport(refBaseAssembly, candidateNonJunctionReads, excludedReads);

            // only add branched assemblies if they have sufficient support
            if(junctionAssembly != assembly)
            {
                if(refBaseAssembly.supportCount() == 0)
                    continue;

                junctionAssembly.mergeRefBaseAssembly(refBaseAssembly);

                if(junctionAssembly.hasUnsetBases())
                    continue;

                // same criteria as above
                if(junctionAssembly.supportCount() < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
                    continue;

                mAssemblies.add(junctionAssembly);
            }
            else
            {
                if(refBaseAssembly.supportCount() > 0) // same criteria as above
                    junctionAssembly.mergeRefBaseAssembly(refBaseAssembly);
            }

            findRemoteRegions(junctionAssembly, excludedReads, discordantReads, remoteJunctionMates, suppJunctionReads);
        }

        // apply filters again since read support may now be below required level

        // set references between them - for now just for TSV output
        for(JunctionAssembly junctionAssembly : mAssemblies)
        {
            mAssemblies.stream().filter(x -> x != junctionAssembly).forEach(x -> junctionAssembly.addBranchedAssembly(x));
        }
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

    private boolean isDiscordantCandidate(final Read read, boolean isForwardJunction, int junctionPosition)
    {
        // cannot cross the junction since will already have considered all junction candidate reads
        if(isForwardJunction)
        {
            if(read.alignmentEnd() >= junctionPosition)
                return false;
        }
        else
        {
            if(read.alignmentStart() <= junctionPosition)
                return false;
        }

        return isDiscordant(read);
    }

    private List<Set<Read>> allocateExcludedReads(final JunctionAssembly assembly, final List<NonJunctionRead> nonJunctionReads)
    {
        // now build out any distinct, branching assemblies from ref-base soft-clips
        List<RefSideSoftClip> refSideSoftClips = assembly.refSideSoftClips();

        boolean isForwardJunction = assembly.isForwardJunction();
        Set<Read> nscReads = Sets.newHashSet();

        int maxRefSideSoftClipPosition = isForwardJunction ?
                refSideSoftClips.stream().filter(x -> !x.matchesOriginal()).mapToInt(x -> x.Position).min().orElse(0)
                : refSideSoftClips.stream().filter(x -> !x.matchesOriginal()).mapToInt(x -> x.Position).max().orElse(0);

        List<Read> allReads = assembly.support().stream().map(x -> x.read()).collect(Collectors.toList());
        nonJunctionReads.forEach(x -> allReads.add(x.read()));

        Set<Read> softClippedReads = Sets.newHashSet();
        refSideSoftClips.stream().filter(x -> !x.matchesOriginal()).forEach(x -> softClippedReads.addAll(x.reads()));

        // find any supporting read or read mate which extends past the furthest soft-clip position
        for(Read read : allReads)
        {
            if(softClippedReads.contains(read) || read.hasMateSet() && softClippedReads.contains(read.mateRead()))
                continue;

            boolean addRead = false;

            if(isForwardJunction)
            {
                if(read.alignmentStart() < maxRefSideSoftClipPosition
                || (read.hasMateSet() && read.mateRead().alignmentStart() < maxRefSideSoftClipPosition))
                {
                    addRead = true;
                }
            }
            else
            {
                if(read.alignmentEnd() > maxRefSideSoftClipPosition
                || (read.hasMateSet() && read.mateRead().alignmentEnd() > maxRefSideSoftClipPosition))
                {
                    addRead = true;
                }
            }

            if(addRead)
            {
                nscReads.add(read);

                if(read.hasMateSet())
                    nscReads.add(read.mateRead());
            }
        }

        int finalAssemblyCount = 1 + refSideSoftClips.size();
        List<Set<Read>> excludedReadsList = Lists.newArrayListWithCapacity(finalAssemblyCount);

        // the original assembly cannot have any
        excludedReadsList.add(softClippedReads);

        for(RefSideSoftClip refSideSoftClip : refSideSoftClips)
        {
            Set<Read> excludedReads = Sets.newHashSet();

            if(!refSideSoftClip.matchesOriginal())
                excludedReads.addAll(nscReads);

            excludedReadsList.add(excludedReads);

            for(RefSideSoftClip other : refSideSoftClips)
            {
                if(refSideSoftClip == other)
                    continue;

                excludedReads.addAll(other.reads());
            }
        }

        return excludedReadsList;
    }

    private static void checkAddRefAssemblySupport(
            final RefBaseAssembly refBaseAssembly, final List<NonJunctionRead> reads, final Set<Read> excludedReads)
    {
        for(NonJunctionRead njRead : reads)
        {
            if(!excludedReads.contains(njRead.read()))
            {
                refBaseAssembly.checkAddRead(
                        njRead.read(), njRead.type(), ASSEMBLY_EXTENSION_BASE_MISMATCH, ASSEMBLY_EXTENSION_OVERLAP_BASES);
            }
        }
    }

    private void findRemoteRegions(
            final JunctionAssembly assembly, final Set<Read> excludedReads,
            final List<Read> discordantReads, final List<Read> remoteJunctionMates, final List<Read> suppJunctionReads)
    {
        if(remoteJunctionMates.isEmpty() && discordantReads.isEmpty())
            return;

        List<RemoteRegion> remoteRegions = Lists.newArrayList();

        discordantReads.stream()
                .filter(x -> !excludedReads.contains(x))
                .forEach(x -> addOrCreateMateRemoteRegion(remoteRegions, x, false));

        for(Read read : remoteJunctionMates)
        {
            if(excludedReads.contains(read))
                continue;

            // look to extend from local mates
            addOrCreateMateRemoteRegion(remoteRegions, read, true);
        }

        for(Read read : suppJunctionReads)
        {
            if(excludedReads.contains(read))
                continue;

            // factor any supplementaries into remote regions
            int scLength = assembly.isForwardJunction() ? read.rightClipLength() : read.leftClipLength();
            addOrCreateSupplementaryRemoteRegion(remoteRegions, read, scLength);
        }

        mergeRegions(remoteRegions);

        // purge regions with only weak supplementary support
        purgeWeakSuppRegions(remoteRegions);

        assembly.addRemoteRegions(remoteRegions);
    }

    private static void addOrCreateMateRemoteRegion(final List<RemoteRegion> remoteRegions, final Read read, boolean isJunctionRead)
    {
        if(read.isMateUnmapped())
            return;

        String mateChr = read.mateChromosome();

        if(mateChr == null || !HumanChromosome.contains(mateChr))
            return;

        addOrCreateRemoteRegion(
                remoteRegions, read, isJunctionRead ? REMOTE_READ_TYPE_JUNCTION_MATE : REMOTE_READ_TYPE_DISCORDANT_READ,
                mateChr, read.mateAlignmentStart(), read.mateAlignmentEnd());
    }

    private static RemoteRegion addOrCreateRemoteRegion(
            final List<RemoteRegion> remoteRegions, final Read read, final int readType,
            final String remoteChr, final int remotePosStart, final int remotePosEnd)
    {
        RemoteRegion matchedRegion = remoteRegions.stream()
                .filter(x -> x.overlaps(remoteChr, remotePosStart, remotePosEnd)).findFirst().orElse(null);

        if(matchedRegion != null)
        {
            matchedRegion.addReadDetails(read.getName(), remotePosStart, remotePosEnd, readType);
            return matchedRegion;
        }
        else
        {
            RemoteRegion newRegion = new RemoteRegion(new ChrBaseRegion(remoteChr, remotePosStart, remotePosEnd), read.getName(), readType);
            remoteRegions.add(newRegion);
            return newRegion;
        }
    }

    private static void addOrCreateSupplementaryRemoteRegion(final List<RemoteRegion> remoteRegions, final Read read, int readScLength)
    {
        SupplementaryReadData suppData = read.supplementaryData();

        if(suppData == null || !HumanChromosome.contains(suppData.Chromosome))
            return;

        int remotePosEnd = getMateAlignmentEnd(suppData.Position, suppData.Cigar);

        /*
        SV_LOGGER.debug("asmJunction({}) read({} flags={}) supp({})",
                mAssembly.junction(), read.getName(), read.getFlags(), suppData);
        */

        RemoteRegion region = addOrCreateRemoteRegion(
                remoteRegions, read, REMOTE_READ_TYPE_JUNCTION_SUPP, suppData.Chromosome, suppData.Position, remotePosEnd);

        region.addSoftClipMapQual(readScLength, suppData.MapQuality);
    }
}
