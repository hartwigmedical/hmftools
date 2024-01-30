package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.getMateAlignmentEnd;
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
import com.hartwig.hmftools.esvee.common.RemoteRegion;
import com.hartwig.hmftools.esvee.read.Read;

public final class AssemblyExtension
{
    public static void extendAssembly(final JunctionAssembly assembly, final List<Read> nonJunctionReads)
    {
        // first establish potential boundaries for extending the assembly on the non-junction side
        List<Read> refExtensionReads = Lists.newArrayList();
        List<Read> discordantReads = Lists.newArrayList();
        List<RemoteRegion> remoteRegions = Lists.newArrayList();

        int minAlignedPosition = assembly.minAlignedPosition();
        int maxAlignedPosition = assembly.maxAlignedPosition();

        boolean isForwardJunction = assembly.initialJunction().isForward();

        // process in order of closest to furtherest out reads in the ref base direction
        List<Read> sortedNonJunctionReads = nonJunctionReads.stream()
                .sorted(Comparator.comparingInt(x -> isForwardJunction ? -x.unclippedEnd() : x.unclippedStart()))
                .collect(Collectors.toList());

        for(Read read : sortedNonJunctionReads)
        {
            if(isForwardJunction)
            {
                // CHECK: must aligned read overlap by any min number of bases to be considered for extension?
                if(read.unclippedEnd() < minAlignedPosition)
                    break;

                minAlignedPosition = min(minAlignedPosition, read.unclippedStart());
            }
            else
            {
                if(read.unclippedStart() > maxAlignedPosition)
                    break;

                maxAlignedPosition = max(maxAlignedPosition, read.unclippedEnd());
            }

            refExtensionReads.add(read);

            // extract any remote positions if this is a discordant read
            if(isDiscordant(read))
            {
                discordantReads.add(read);

                addOrCreateMateRemoteRegion(remoteRegions, read);
            }

            // alternatively if it is local it may have a mate even further away, but is this interesting?
        }

        // add remote regions from junction read supplementaries
        for(AssemblySupport support : assembly.support())
        {
            // look to extend from local mates
            if(isDiscordant(support.read()))
            {
                // discordantReads.add(support.read());

                addOrCreateMateRemoteRegion(remoteRegions, support.read());
            }

            // factor any supplementaries into remote regions
            addOrCreateSupplementaryRemoteRegion(remoteRegions, support.read());
        }

        RemoteRegion.mergeRegions(remoteRegions);

        assembly.addRemoteRegions(remoteRegions);
    }

    private static void addOrCreateMateRemoteRegion(final List<RemoteRegion> remoteRegions, final Read read)
    {
        if(read.isMateUnmapped())
            return;

        String mateChr = read.mateChromosome();

        if(!HumanChromosome.contains(mateChr))
            return;

        addOrCreateRemoteRegion(remoteRegions, read, mateChr, read.mateAlignmentStart(), read.mateAlignmentEnd());
    }

    private static void addOrCreateRemoteRegion(
            final List<RemoteRegion> remoteRegions, final Read read, final String remoteChr, final int remotePosStart, final int remotePosEnd)
    {
        RemoteRegion matchedRegion = remoteRegions.stream()
                .filter(x -> x.overlaps(remoteChr, remotePosStart, remotePosEnd)).findFirst().orElse(null);

        if(matchedRegion != null)
        {
            matchedRegion.addReadDetails(read.getName(), remotePosStart, remotePosEnd);
        }
        else
        {
            remoteRegions.add(new RemoteRegion(new ChrBaseRegion(remoteChr, remotePosStart, remotePosEnd), read.getName()));
        }

    }

    private static void addOrCreateSupplementaryRemoteRegion(final List<RemoteRegion> remoteRegions, final Read read)
    {
        SupplementaryReadData suppData = read.supplementaryData();

        if(suppData == null || !HumanChromosome.contains(suppData.Chromosome))
            return;

        int remotePosEnd = getMateAlignmentEnd(suppData.Position, suppData.Cigar);

        addOrCreateRemoteRegion(remoteRegions, read, suppData.Chromosome, suppData.Position, remotePosEnd);
    }
}
