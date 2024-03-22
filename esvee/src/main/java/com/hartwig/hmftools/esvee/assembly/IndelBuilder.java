package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.esvee.AssemblyConstants.MIN_INDEL_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.RefBaseExtender.isValidSupportCoordsVsJunction;
import static com.hartwig.hmftools.esvee.types.SupportType.JUNCTION_MATE;
import static com.hartwig.hmftools.esvee.read.ReadUtils.isDiscordant;

import static htsjdk.samtools.CigarOperator.I;

import java.util.List;

import com.hartwig.hmftools.esvee.types.AssemblySupport;
import com.hartwig.hmftools.esvee.types.IndelCoords;
import com.hartwig.hmftools.esvee.types.Junction;
import com.hartwig.hmftools.esvee.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.read.ReadFilters;

import htsjdk.samtools.CigarElement;

public final class IndelBuilder
{
    public static void findIndelExtensionReads(
            final Junction junction, final List<Read> rawReads,
            final List<Read> extensionReads, final List<Read> junctionReads, final List<Read> nonJunctionReads)
    {
        for(Read read : rawReads)
        {
            IndelCoords indelCoords = read.indelCoords();

            if(indelCoords != null)
            {
                // must match junction exactly to be considered for support
                if(!indelCoords.matchesJunction(junction.Position, junction.Orientation))
                    continue;

                if(indelCoords.Length >= MIN_INDEL_LENGTH)
                    extensionReads.add(read);
                else
                    junctionReads.add(read);
            }
            else
            {
                if(!ReadFilters.recordSoftClipsAndCrossesJunction(read, junction))
                {
                    nonJunctionReads.add(read);
                    continue;
                }

                junctionReads.add(read);
            }
        }
    }

    public static boolean convertedIndelCrossesJunction(final JunctionAssembly assembly, final Read read)
    {
        if(assembly.isForwardJunction())
        {
            if(read.unclippedEnd() > assembly.junction().Position)
                return read.indelImpliedAlignmentEnd() > 0 || read.isConvertedIndel();
        }
        else
        {
            if(read.unclippedStart() < assembly.junction().Position)
                return read.indelImpliedAlignmentStart() > 0 || read.isConvertedIndel();
        }

        return false;
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

            // no discordant candidates for local indels
            if(isJunctionMate)
            {
                assembly.addCandidateSupport(read, JUNCTION_MATE);
            }
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
