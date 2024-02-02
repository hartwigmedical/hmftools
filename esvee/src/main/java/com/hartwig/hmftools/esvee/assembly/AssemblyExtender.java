package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.getMateAlignmentEnd;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.ASSEMBLY_EXTENSION_BASE_MISMATCH;
import static com.hartwig.hmftools.esvee.SvConstants.ASSEMBLY_EXTENSION_OVERLAP_BASES;
import static com.hartwig.hmftools.esvee.common.RemoteRegion.REMOTE_READ_TYPE_DISCORDANT_READ;
import static com.hartwig.hmftools.esvee.common.RemoteRegion.REMOTE_READ_TYPE_JUNCTION_MATE;
import static com.hartwig.hmftools.esvee.common.RemoteRegion.REMOTE_READ_TYPE_JUNCTION_SUPP;
import static com.hartwig.hmftools.esvee.common.SupportType.DISCORDANT;
import static com.hartwig.hmftools.esvee.common.SupportType.JUNCTION_MATE;
import static com.hartwig.hmftools.esvee.read.ReadUtils.isDiscordant;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.esvee.common.AssemblySupport;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.RefBaseAssembly;
import com.hartwig.hmftools.esvee.common.RemoteRegion;
import com.hartwig.hmftools.esvee.common.SupportType;
import com.hartwig.hmftools.esvee.read.Read;

public class AssemblyExtender
{
    private final JunctionAssembly mAssembly;

    public AssemblyExtender(final JunctionAssembly assembly)
    {
        mAssembly = assembly;
    }

    public void extendAssembly(final List<Read> nonJunctionReads)
    {
        // first establish potential boundaries for extending the assembly on the non-junction side
        int minAlignedPosition = mAssembly.minAlignedPosition();
        int maxAlignedPosition = mAssembly.maxAlignedPosition();

        boolean isForwardJunction = mAssembly.junction().isForward();
        int junctionPosition = mAssembly.junction().Position;

        // process in order of closest to furthest-out reads in the ref base direction
        List<Read> discordantReads = nonJunctionReads.stream()
                .filter(x -> isDiscordant(x) || x.isMateUnmapped())
                .filter(x -> !x.hasMateSet() || !mAssembly.hasReadSupport(x.mateRead()))
                .collect(Collectors.toList());

        List<NonJunctionRead> candidateReads = discordantReads.stream().map(x -> new NonJunctionRead(x, DISCORDANT)).collect(Collectors.toList());

        List<Read> remoteJunctionMates = Lists.newArrayList();
        List<Read> suppJunctionReads = Lists.newArrayList();

        // add any junction mates in the same window
        for(AssemblySupport support : mAssembly.support())
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

            if(mateRead == null || discordantReads.contains(mateRead))
                continue;

            if(isForwardJunction)
            {
                if(mateRead.unclippedEnd() >= junctionPosition)
                    continue;
            }
            else
            {
                if(mateRead.unclippedStart() <= junctionPosition)
                    continue;
            }

            candidateReads.add(new NonJunctionRead(mateRead, JUNCTION_MATE));
        }

        findRemoteRegions(discordantReads, remoteJunctionMates, suppJunctionReads);

        if(candidateReads.isEmpty())
            return;

        List<NonJunctionRead> sortedNonJunctionReads = candidateReads.stream()
                .sorted(Comparator.comparingInt(x -> isForwardJunction ? -x.read().unclippedEnd() : x.read().unclippedStart()))
                .collect(Collectors.toList());

        for(NonJunctionRead njRead : sortedNonJunctionReads)
        {
            Read read = njRead.read();

            if(isForwardJunction)
            {
                if(read.unclippedEnd() < minAlignedPosition + ASSEMBLY_EXTENSION_OVERLAP_BASES)
                    break;

                minAlignedPosition = min(minAlignedPosition, read.unclippedStart());
            }
            else
            {
                if(read.unclippedStart() > maxAlignedPosition - ASSEMBLY_EXTENSION_OVERLAP_BASES)
                    break;

                maxAlignedPosition = max(maxAlignedPosition, read.unclippedEnd());
            }
        }

        // find a support read which extend out the furthest from the junction and with the least mismatches
        int extensionRefPosition = isForwardJunction ? minAlignedPosition : maxAlignedPosition;
        RefBaseAssembly refBaseAssembly = new RefBaseAssembly(mAssembly, extensionRefPosition);

        for(NonJunctionRead njRead : sortedNonJunctionReads)
        {
            refBaseAssembly.checkAddRead(njRead.read(), njRead.type(), ASSEMBLY_EXTENSION_BASE_MISMATCH);
        }

        mAssembly.setRefBaseAssembly(refBaseAssembly);
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
    }

    private void findRemoteRegions(final List<Read> discordantReads, final List<Read> remoteJunctionMates, final List<Read> suppJunctionReads)
    {
        if(remoteJunctionMates.isEmpty() && discordantReads.isEmpty())
            return;

        List<RemoteRegion> remoteRegions = Lists.newArrayList();

        discordantReads.forEach(x -> addOrCreateMateRemoteRegion(remoteRegions, x, false));

        for(Read read : remoteJunctionMates)
        {
            // look to extend from local mates
            addOrCreateMateRemoteRegion(remoteRegions, read, true);
        }

        for(Read read : suppJunctionReads)
        {
            // factor any supplementaries into remote regions
            addOrCreateSupplementaryRemoteRegion(remoteRegions, read);
        }

        RemoteRegion.mergeRegions(remoteRegions);

        mAssembly.addRemoteRegions(remoteRegions);
    }

    private static void addOrCreateMateRemoteRegion(final List<RemoteRegion> remoteRegions, final Read read, boolean isJunctionRead)
    {
        if(read.isMateUnmapped())
            return;

        String mateChr = read.mateChromosome();

        if(!HumanChromosome.contains(mateChr))
            return;

        addOrCreateRemoteRegion(
                remoteRegions, read, isJunctionRead ? REMOTE_READ_TYPE_JUNCTION_MATE : REMOTE_READ_TYPE_DISCORDANT_READ,
                mateChr, read.mateAlignmentStart(), read.mateAlignmentEnd());
    }

    private static void addOrCreateRemoteRegion(
            final List<RemoteRegion> remoteRegions, final Read read, final int readType,
            final String remoteChr, final int remotePosStart, final int remotePosEnd)
    {
        RemoteRegion matchedRegion = remoteRegions.stream()
                .filter(x -> x.overlaps(remoteChr, remotePosStart, remotePosEnd)).findFirst().orElse(null);

        if(matchedRegion != null)
        {
            matchedRegion.addReadDetails(read.getName(), remotePosStart, remotePosEnd, readType);
        }
        else
        {
            remoteRegions.add(new RemoteRegion(new ChrBaseRegion(remoteChr, remotePosStart, remotePosEnd), read.getName(), readType));
        }
    }

    private void addOrCreateSupplementaryRemoteRegion(final List<RemoteRegion> remoteRegions, final Read read)
    {
        SupplementaryReadData suppData = read.supplementaryData();

        if(suppData == null || !HumanChromosome.contains(suppData.Chromosome))
            return;

        int remotePosEnd = getMateAlignmentEnd(suppData.Position, suppData.Cigar);

        /*
        SV_LOGGER.debug("asmJunction({}) read({} flags={}) supp({})",
                mAssembly.junction(), read.getName(), read.getFlags(), suppData);
        */

        addOrCreateRemoteRegion(remoteRegions, read, REMOTE_READ_TYPE_JUNCTION_SUPP, suppData.Chromosome, suppData.Position, remotePosEnd);
    }
}
