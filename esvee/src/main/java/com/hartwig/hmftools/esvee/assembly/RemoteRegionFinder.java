package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.getMateAlignmentEnd;
import static com.hartwig.hmftools.esvee.common.RemoteRegion.REMOTE_READ_TYPE_DISCORDANT_READ;
import static com.hartwig.hmftools.esvee.common.RemoteRegion.REMOTE_READ_TYPE_JUNCTION_MATE;
import static com.hartwig.hmftools.esvee.common.RemoteRegion.REMOTE_READ_TYPE_JUNCTION_SUPP;
import static com.hartwig.hmftools.esvee.common.RemoteRegion.mergeRegions;
import static com.hartwig.hmftools.esvee.common.RemoteRegion.purgeWeakSupplementaryRegions;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.RemoteRegion;
import com.hartwig.hmftools.esvee.read.Read;

public final class RemoteRegionFinder
{

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
