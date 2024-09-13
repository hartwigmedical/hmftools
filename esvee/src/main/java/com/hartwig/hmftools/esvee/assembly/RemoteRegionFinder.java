package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.getMateAlignmentEnd;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.types.RemoteReadType.DISCORDANT;
import static com.hartwig.hmftools.esvee.assembly.types.RemoteReadType.JUNCTION_MATE;
import static com.hartwig.hmftools.esvee.assembly.types.RemoteReadType.SUPPLEMENTARY;
import static com.hartwig.hmftools.esvee.assembly.types.RemoteRegion.mergeRegions;
import static com.hartwig.hmftools.esvee.assembly.types.RemoteRegion.purgeWeakSupplementaryRegions;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.RemoteReadType;
import com.hartwig.hmftools.esvee.assembly.types.RemoteRegion;
import com.hartwig.hmftools.esvee.assembly.read.Read;

public final class RemoteRegionFinder
{
    private static final int MAX_REMOTE_REGION_LENGTH = 10000;

    public static void findRemoteRegions(
            final JunctionAssembly assembly, final List<Read> discordantReads, final List<Read> remoteJunctionMates,
            final List<Read> suppJunctionReads)
    {
        if(remoteJunctionMates.isEmpty() && discordantReads.isEmpty() && suppJunctionReads.isEmpty())
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
            int scLength = assembly.isForwardJunction() ? read.rightClipLength() : read.leftClipLength();
            addOrCreateSupplementaryRemoteRegion(remoteRegions, read, scLength);
        }

        mergeRegions(remoteRegions);

        // purge regions with only weak supplementary support
        purgeWeakSupplementaryRegions(remoteRegions);

        // remove regions which overlap the assembly
        int index = 0;
        while(index < remoteRegions.size())
        {
            RemoteRegion remoteRegion = remoteRegions.get(index);

            if(assemblyOverlapsRemoteRegion(assembly, remoteRegion))
            {
                remoteRegions.remove(index);
            }
            else if(remoteRegion.length() > MAX_REMOTE_REGION_LENGTH)
            {
                SV_LOGGER.warn("assembly({}) excluding long remote region({})", assembly, remoteRegion);
                remoteRegions.remove(index);
            }
            else
            {
                ++index;
            }
        }

        assembly.addRemoteRegions(remoteRegions);
    }

    private static boolean assemblyOverlapsRemoteRegion(final JunctionAssembly assembly, final RemoteRegion remoteRegion)
    {
        if(!assembly.junction().Chromosome.equals(remoteRegion.Chromosome))
            return false;

        if(assembly.isForwardJunction())
            return positionsOverlap(assembly.refBasePosition(), assembly.junction().Position, remoteRegion.start(), remoteRegion.end());
        else
            return positionsOverlap(assembly.junction().Position, assembly.refBasePosition(), remoteRegion.start(), remoteRegion.end());
    }

    private static void addOrCreateMateRemoteRegion(final List<RemoteRegion> remoteRegions, final Read read, boolean isJunctionRead)
    {
        if(read.isMateUnmapped() || read.isSupplementary())
            return;

        String mateChr = read.mateChromosome();

        if(mateChr == null || !HumanChromosome.contains(mateChr))
            return;

        addOrCreateRemoteRegion(
                remoteRegions, read, isJunctionRead ? JUNCTION_MATE : DISCORDANT,
                mateChr, read.mateAlignmentStart(), read.mateAlignmentEnd());
    }

    private static RemoteRegion addOrCreateRemoteRegion(
            final List<RemoteRegion> remoteRegions, final Read read, final RemoteReadType readType,
            final String remoteChr, final int remotePosStart, final int remotePosEnd)
    {
        RemoteRegion matchedRegion = remoteRegions.stream()
                .filter(x -> x.overlaps(remoteChr, remotePosStart, remotePosEnd)).findFirst().orElse(null);

        if(matchedRegion != null)
        {
            matchedRegion.addReadDetails(read.id(), remotePosStart, remotePosEnd, readType);
            return matchedRegion;
        }
        else
        {
            RemoteRegion newRegion = new RemoteRegion(
                    new ChrBaseRegion(remoteChr, remotePosStart, remotePosEnd), read.id(), readType);
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

        RemoteRegion region = addOrCreateRemoteRegion(
                remoteRegions, read, SUPPLEMENTARY, suppData.Chromosome, suppData.Position, remotePosEnd);

        region.addSoftClipMapQual(readScLength, suppData.MapQuality);
    }
}
