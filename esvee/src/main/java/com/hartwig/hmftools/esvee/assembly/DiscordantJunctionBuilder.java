package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_DISCORDANT_MIN_MAP_QUALITY;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.MAX_OBSERVED_CONCORDANT_FRAG_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.RemoteRegionFinder.addOrCreateMateRemoteRegion;
import static com.hartwig.hmftools.esvee.assembly.types.RemoteRegion.mergeRegions;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.RemoteRegion;

public class DiscordantJunctionBuilder
{
    private final Junction mJunction;
    private Junction mNewJunction;
    private final List<Read> mNonJunctionReads;

    public DiscordantJunctionBuilder(final Junction junction, final List<Read> nonJunctionReads)
    {
        mJunction = junction;
        mNewJunction = null;
        mNonJunctionReads = nonJunctionReads;
    }

    public boolean hasNewJunction() { return mNewJunction != null; }
    public Junction newJunction() { return mNewJunction; }

    public void assessDiscordantJunction(final List<Read> rawReads, final List<Read> extensionReads, final List<Read> junctionReads)
    {
        if(rawReads.size() < ASSEMBLY_MIN_READ_SUPPORT)
            return;

        int originalJuncPosition = mJunction.Position;

        // find the applicable remote region from amongst
        List<Read> candidateJunctionReads = rawReads.stream()
                .filter(x -> (mJunction.isForward() && x.alignmentEnd() == originalJuncPosition)
                        || (mJunction.isReverse() && x.alignmentStart() == originalJuncPosition))
                .collect(Collectors.toList());

        if(candidateJunctionReads.isEmpty())
            return;

        List<RemoteRegion> remoteRegions = Lists.newArrayList();
        candidateJunctionReads.forEach(x -> addOrCreateMateRemoteRegion(remoteRegions, x, true));

        if(remoteRegions.isEmpty())
            return;

        if(remoteRegions.size() > 1)
        {
            mergeRegions(remoteRegions);
            Collections.sort(remoteRegions, Comparator.comparingInt(x -> -x.readCount()));
        }

        RemoteRegion mainRemoteRegion = remoteRegions.get(0);

        // take the inner-most read position and its remote region, then only use other reads with matching remote regions
        List<Integer> extensionJuncPositions = Lists.newArrayListWithCapacity(rawReads.size());

        Set<String> candidateReadIds = Sets.newHashSet();

        int minReadPosStart, maxReadPosEnd;

        if(mJunction.isForward())
        {
            minReadPosStart = max(1, originalJuncPosition - MAX_OBSERVED_CONCORDANT_FRAG_LENGTH);
            maxReadPosEnd = originalJuncPosition;
        }
        else
        {
            minReadPosStart = originalJuncPosition;
            maxReadPosEnd = originalJuncPosition + MAX_OBSERVED_CONCORDANT_FRAG_LENGTH;
        }

        for(Read read : rawReads)
        {
            if(read.mappingQuality() < ASSEMBLY_DISCORDANT_MIN_MAP_QUALITY)
                continue;

            // check that the read maps to the same remote region to qualify as an junction candidate
            if(!mainRemoteRegion.overlaps(read.mateChromosome(), read.mateAlignmentStart(), read.mateAlignmentEnd()))
                continue;

            // ensure reads aren't past the inner-most discordant read, nor outside valid bounds
            if(read.alignmentStart() < minReadPosStart || read.alignmentEnd() > maxReadPosEnd)
                continue;

            candidateReadIds.add(read.id());

            if(mJunction.isForward())
            {
                if(read.alignmentEnd() > originalJuncPosition)
                    continue;

                extensionJuncPositions.add(read.unclippedEnd());
            }
            else
            {
                if(read.alignmentStart() < originalJuncPosition)
                    continue;

                extensionJuncPositions.add(read.unclippedStart());
            }
        }

        if(extensionJuncPositions.size() < ASSEMBLY_MIN_READ_SUPPORT)
            return;

        if(mJunction.isForward())
            Collections.sort(extensionJuncPositions, Comparator.reverseOrder());
        else
            Collections.sort(extensionJuncPositions);

        int outerJuncPosition = extensionJuncPositions.get(0);
        int secondJuncPosition = extensionJuncPositions.get(ASSEMBLY_MIN_READ_SUPPORT - 1);
        int adjustedJuncPosition;

        if(mJunction.isForward())
        {
            outerJuncPosition -= ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
            secondJuncPosition -= ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;
            adjustedJuncPosition = min(outerJuncPosition, secondJuncPosition);
        }
        else
        {
            outerJuncPosition += ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
            secondJuncPosition += ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;
            adjustedJuncPosition = max(outerJuncPosition, secondJuncPosition);
        }

        mNewJunction = new Junction(
                mJunction.Chromosome, adjustedJuncPosition, mJunction.Orient, true, false, false,
                null); // no SAGA-match owing to imprecise junction coords

        mNewJunction.setRawDiscordantPosition(originalJuncPosition);

        for(Read read : rawReads)
        {
            if(read.alignmentStart() < minReadPosStart || read.alignmentEnd() > maxReadPosEnd)
                continue;

            // skip indels supporting discordant-only junctions
            if(mJunction.isForward())
            {
                if(read.hasIndelImpliedUnclippedEnd())
                    continue;
            }
            else
            {
                if(read.hasIndelImpliedUnclippedStart())
                    continue;
            }

            if(!candidateReadIds.contains(read.id()) || read.mappingQuality() < ASSEMBLY_DISCORDANT_MIN_MAP_QUALITY)
            {
                mNonJunctionReads.add(read);
                continue;
            }

            int extensionLength;

            if(mJunction.isForward())
            {
                if(read.alignmentEnd() > originalJuncPosition)
                    continue;

                extensionLength = read.alignmentEnd() >= adjustedJuncPosition ? read.unclippedEnd() - adjustedJuncPosition : -1;
            }
            else
            {
                if(read.alignmentStart() < originalJuncPosition)
                    continue;

                extensionLength = read.alignmentStart() <= adjustedJuncPosition ? adjustedJuncPosition - read.unclippedStart() : -1;
            }

            if(extensionLength >= ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH)
                extensionReads.add(read);
            else if(extensionLength > 0)
                junctionReads.add(read);
            else
                mNonJunctionReads.add(read);
        }
    }
}
