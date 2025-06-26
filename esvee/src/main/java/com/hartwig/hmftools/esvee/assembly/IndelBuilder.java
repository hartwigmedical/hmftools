package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_INDEL_UNLINKED_ASSEMBLY_INDEL_PERC;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_INDEL_UNLINKED_ASSEMBLY_MIN_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.INDEL_TO_SC_MIN_SIZE_SOFTCLIP;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.NO_LINK;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.UNSET;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_LENGTH;

import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;
import com.hartwig.hmftools.esvee.common.IndelCoords;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.read.Read;

import htsjdk.samtools.CigarElement;

public final class IndelBuilder
{
    public static boolean hasIndelJunctionReads(final Junction junction, final List<Read> reads)
    {
        return reads.stream()
                .filter(x -> x.indelCoords() != null && x.indelCoords().Length >= MIN_INDEL_LENGTH)
                .anyMatch(x -> x.indelCoords().matchesJunction(junction.Position, junction.Orient));
    }

    public static void findIndelExtensionReads(
            final Junction junction, final List<Read> rawReads,
            final List<Read> extensionReads, final List<Read> junctionReads, final List<Read> nonJunctionReads)
    {
        Map<Integer,Integer> indelCoordsCount = Maps.newHashMap();
        IndelCoords dominantIndel = null;
        int dominantIndelFreq = 0;

        for(Read read : rawReads)
        {
            if(read.invalidIndel())
                continue;

            IndelCoords indelCoords = read.indelCoords();

            if(indelCoords != null)
            {
                // must match junction exactly on that side to be considered for support
                if(!indelCoords.matchesJunction(junction.Position, junction.Orient))
                    continue;

                if(indelCoords.Length >= MIN_INDEL_LENGTH)
                {
                    extensionReads.add(read);

                    int indelLengthFreq = indelCoordsCount.getOrDefault(indelCoords.Length, 0) + 1;
                    indelCoordsCount.put(indelCoords.Length, indelLengthFreq);

                    if(indelLengthFreq > dominantIndelFreq)
                    {
                        dominantIndel = indelCoords;
                        dominantIndelFreq = indelLengthFreq;
                    }
                }
                else
                {
                    junctionReads.add(read);
                }
            }
            else
            {
                // testing shorter indel and soft-clipped reads
                boolean crossesJunction = false;

                if(junction.isForward())
                {
                    if(read.hasIndelImpliedUnclippedEnd() && read.maxUnclippedEnd() > junction.Position)
                    {
                        crossesJunction = true;
                    }
                    else if(read.isRightClipped() && read.unclippedEnd() > junction.Position)
                    {
                        crossesJunction = true;
                    }
                }
                else
                {
                    if(read.hasIndelImpliedUnclippedStart() && read.minUnclippedStart() < junction.Position)
                    {
                        crossesJunction = true;
                    }
                    else if(read.isLeftClipped() && read.unclippedStart() < junction.Position)
                    {
                        crossesJunction = true;
                    }
                }

                if(crossesJunction)
                    junctionReads.add(read);
                else
                    nonJunctionReads.add(read);
            }
        }

        if(dominantIndel != null)
            junction.setIndelCoords(dominantIndel);
    }

    public static boolean calcIndelInferredUnclippedPositions(final Read read)
    {
        if(read.cigarElements().size() < 3)
            return false;

        int insertTotal = 0;
        int alignedTotal = 0;

        boolean hasMinLengthIndel = false;
        int lastIndex = read.cigarElements().size() - 1;

        for(int i = 0; i < read.cigarElements().size(); ++i)
        {
            CigarElement element = read.cigarElements().get(i);

            if(element.getOperator().isIndel())
            {
                hasMinLengthIndel |= element.getLength() >= INDEL_TO_SC_MIN_SIZE_SOFTCLIP;

                // must be straddled by aligned bases
                if(i == 0 || (i == 1 && read.cigarElements().get(0).getOperator() != M))
                {
                    read.markInvalidIndel();
                    return false;
                }

                if(i == lastIndex || (i == lastIndex - 1 && read.cigarElements().get(lastIndex).getOperator() != M))
                {
                    read.markInvalidIndel();
                    return false;
                }

                if(element.getOperator() == I)
                    insertTotal += element.getLength();
            }
            else if(element.getOperator() == M)
            {
                alignedTotal += element.getLength();
            }
        }

        if(!hasMinLengthIndel)
            return false;

        Integer unclippedStart = null;
        Integer unclippedEnd = null;

        if(!read.isLeftClipped())
            unclippedStart = read.alignmentEnd() - alignedTotal + 1 - insertTotal;

        if(!read.isRightClipped())
            unclippedEnd = read.alignmentStart() + alignedTotal - 1 + insertTotal;

        read.setIndelInferredUnclippedPositions(unclippedStart, unclippedEnd);
        return true;
    }

    public static boolean convertedIndelCrossesJunction(final Junction junction, final Read read)
    {
        if(junction.isForward())
        {
            if(read.unclippedEnd() > junction.Position)
                return read.hasIndelImpliedUnclippedEnd();
        }
        else
        {
            if(read.unclippedStart() < junction.Position)
                return read.hasIndelImpliedUnclippedStart();
        }

        return false;
    }

    public static boolean indelAssembliesMatch(final JunctionAssembly first, final JunctionAssembly second)
    {
        if(first.indelCoords() == null || second.indelCoords() == null || first.isForwardJunction() == second.isForwardJunction())
            return false;

        return first.indelCoords().matches(second.indelCoords());
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

    public static boolean isWeakIndelBasedUnlinkedAssembly(final JunctionAssembly assembly)
    {
        if(assembly.junction().indelBased()) // only applicable to soft-clipped assemblies
            return false;

        if(assembly.outcome() != UNSET && assembly.outcome() != NO_LINK)
            return false;

        if(assembly.extensionLength() >= ASSEMBLY_INDEL_UNLINKED_ASSEMBLY_MIN_LENGTH)
            return false;

        // classify as weak if has the top 2 indel reads have the longest extensions
        int totalJuncReads = 0;
        int maxNonIndelLength = 0;
        List<Integer> indelReadLengths = Lists.newArrayList();

        for(SupportRead read : assembly.support())
        {
            if(read.type().isSplitSupport())
            {
                int extensionLength = read.extensionBaseMatches();
                ++totalJuncReads;

                if(read.type() == SupportType.INDEL)
                {
                    indelReadLengths.add(extensionLength);
                }
                else
                {
                    maxNonIndelLength = max(maxNonIndelLength, extensionLength);
                }
            }
        }

        if(indelReadLengths.size() >= ASSEMBLY_INDEL_UNLINKED_ASSEMBLY_INDEL_PERC * totalJuncReads)
            return true;

        if(indelReadLengths.size() < 2)
            return false;

        Collections.sort(indelReadLengths);
        int secondLongestIndel = indelReadLengths.get(indelReadLengths.size() - 2);

        return secondLongestIndel > maxNonIndelLength;
    }
}
