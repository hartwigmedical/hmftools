package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MAX_BASE_MISMATCH;
import static com.hartwig.hmftools.esvee.assembly.RefBaseExtender.isValidSupportCoordsVsJunction;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.readQualFromJunction;
import static com.hartwig.hmftools.esvee.common.SupportType.CANDIDATE_DISCORDANT;
import static com.hartwig.hmftools.esvee.common.SupportType.INDEL;
import static com.hartwig.hmftools.esvee.common.SupportType.JUNCTION_MATE;
import static com.hartwig.hmftools.esvee.read.ReadUtils.isDiscordant;

import static htsjdk.samtools.CigarOperator.I;

import java.util.List;
import java.util.Map;
import java.util.Optional;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.common.AssemblySupport;
import com.hartwig.hmftools.esvee.common.IndelCoords;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.read.Read;

import htsjdk.samtools.CigarElement;

public final class IndelBuilder
{
    protected static JunctionAssembly buildFromIndelReads(
            final Junction junction, final List<Read> indelReads, final List<Read> shortIndelReads, final List<Read> softClippedReads)
    {
        // find the most common indel length and the read with the highest qual bases from the junction
        Map<Integer,Integer> lengthFrequency = Maps.newHashMap();

        for(Read read : indelReads)
        {
            IndelCoords indelCoords  = read.indelCoords();
            lengthFrequency.compute(indelCoords.Length, (key, oldValue) -> Optional.ofNullable(oldValue).orElse(0) + 1);
        }

        // find the longest indel with the most frequency
        int maxIndelLength = 0;
        int maxIndelLengthFreq = 0;

        for(Map.Entry<Integer,Integer> entry : lengthFrequency.entrySet())
        {
            if(entry.getValue() > maxIndelLengthFreq || (entry.getValue() == maxIndelLengthFreq && maxIndelLength > entry.getKey()))
            {
                maxIndelLengthFreq = entry.getValue();
                maxIndelLength = entry.getKey();
            }
        }

        Read maxJunctionBaseQualRead = null;
        int maxJunctionBaseQualTotal = 0;

        // now the base indel coords have been selected, look for the max extension only amongst reads that match that
        int minAlignedPosition = junction.Position;
        int maxAlignedPosition = junction.Position;
        int maxDistanceFromJunction = 0;

        for(Read read : indelReads)
        {
            if(read.indelCoords().Length != maxIndelLength)
                continue;

            int readJunctionIndex = read.getReadIndexAtReferencePosition(junction.Position, true);

            int extensionDistance = junction.isForward() ? read.basesLength() - readJunctionIndex - 1 : readJunctionIndex;

            maxDistanceFromJunction = max(maxDistanceFromJunction, extensionDistance);

            if(junction.isForward())
            {
                maxAlignedPosition = max(maxAlignedPosition, read.unclippedEnd());
            }
            else
            {
                minAlignedPosition = min(minAlignedPosition, read.unclippedStart());
            }

            int junctionBaseQualTotal = readQualFromJunction(read, junction);

            if(junctionBaseQualTotal > maxJunctionBaseQualTotal)
            {
                maxJunctionBaseQualTotal = junctionBaseQualTotal;
                maxJunctionBaseQualRead = read;
            }
        }

        JunctionAssembly assembly = new JunctionAssembly(
                junction, maxJunctionBaseQualRead, maxDistanceFromJunction, minAlignedPosition,  maxAlignedPosition);

        // assembly.markIndel(); // uses junction attribute from SvPrep

        int permittedMismatches = PRIMARY_ASSEMBLY_MAX_BASE_MISMATCH;

        for(Read read : indelReads)
        {
            if(read == maxJunctionBaseQualRead)
                continue;

            if(read.indelCoords().Length != maxIndelLength) //  || !assembly.checkReadMatches(read, permittedMismatches)
                continue;

            assembly.addReadSupport(read, INDEL, true);
        }

        for(Read read : shortIndelReads)
        {
            if(!assembly.checkReadMatches(read, permittedMismatches))
                continue;

            assembly.addJunctionRead(read, true);
        }

        for(Read read : softClippedReads)
        {
            if(!assembly.checkReadMatches(read, permittedMismatches))
                continue;

            assembly.addJunctionRead(read, true);
        }

        return assembly;
    }

    public static void findIndelExtensions(final JunctionAssembly assembly, final List<Read> unfilteredNonJunctionReads)
    {
        // add junction mates only, could consider add reads which span since these should have a corresponding read in the other junction
        final IndelCoords indelCoords = assembly.initialRead().indelCoords();
        boolean isForwardJunction = assembly.junction().isForward();
        int junctionPosition = assembly.junction().Position;

        int indelInnerStart = indelCoords.PosStart + 1;
        int indelInnerEnd = indelCoords.PosEnd - 1;

        for(Read read : unfilteredNonJunctionReads)
        {
            boolean isJunctionMate = false;
            boolean alreadySupport = false;

            // check vs side of the junction and its orientation
            if(!isValidSupportCoordsVsJunction(read, isForwardJunction, junctionPosition))
                continue;

            for(AssemblySupport support : assembly.support())
            {
                if(support.read() == read)
                {
                    alreadySupport = true;
                    break;
                }

                if(support.read().mateRead() == read)
                {
                    isJunctionMate = true;
                    break;
                }
            }

            if(alreadySupport)
                continue;

            if(positionsOverlap(indelInnerStart, indelInnerEnd, read.alignmentStart(), read.alignmentEnd()))
            {
                continue;
            }

            if(!isJunctionMate)
            {
                if(isDiscordant(read) || read.isMateUnmapped())
                    continue;

                if(positionsOverlap(indelInnerStart, indelInnerEnd, read.mateAlignmentStart(), read.mateAlignmentEnd()))
                {
                    continue;
                }

                // 'discordant' reads must have their mate on the other side of the indel
                boolean onLowerSide = read.alignmentEnd() <= indelCoords.PosStart;
                boolean mateOnLowerSide = read.mateAlignmentEnd() <= indelCoords.PosStart;

                if(onLowerSide == mateOnLowerSide)
                    continue;
            }

            assembly.addCandidateSupport(read, isJunctionMate ? JUNCTION_MATE : CANDIDATE_DISCORDANT);
        }
    }

    public static String findInsertedBases(final Read read)
    {
        final IndelCoords indelCoords = read.indelCoords();

        int readIndex = 0;
        int refPosition = read.alignmentStart();

        for(CigarElement element : read.cigarElements())
        {
            if(element.getOperator() == I && element.getLength() == indelCoords.Length
            || refPosition == indelCoords.PosEnd)
            {
                StringBuilder insertedBases = new StringBuilder();

                for(int i = 0; i < element.getLength(); ++i)
                {
                    insertedBases.append((char)read.getBases()[readIndex + i]);
                }

                return insertedBases.toString();
            }

            if(element.getOperator().consumesReadBases())
                readIndex += element.getLength();

            if(element.getOperator().consumesReferenceBases())
                refPosition += element.getLength();
        }

        return "";
    }
}
