package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.bam.CigarUtils.getPositionFromReadIndex;
import static com.hartwig.hmftools.common.bam.CigarUtils.maxIndelLength;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.INDEL_TO_SC_MAX_SIZE_SOFTCLIP;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.INDEL_TO_SC_MIN_SIZE_SOFTCLIP;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_LENGTH;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;
import com.hartwig.hmftools.esvee.common.IndelCoords;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.read.ReadFilters;

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
        for(Read read : rawReads)
        {
            IndelCoords indelCoords = read.indelCoords();

            if(indelCoords != null)
            {
                // must match junction exactly to be considered for support
                if(!indelCoords.matchesJunction(junction.Position, junction.Orient))
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

    public static boolean calcIndelInferredUnclippedPositions(final Read read)
    {
        if(read.cigarElements().size() < 3)
            return false;

        int insertTotal = 0;
        int alignedTotal = 0;

        boolean hasMinLengthIndel = false;

        for(int i = 0; i < read.cigarElements().size(); ++i)
        {
            CigarElement element = read.cigarElements().get(i);

            switch(element.getOperator())
            {
                case D:
                    hasMinLengthIndel |= element.getLength() >= INDEL_TO_SC_MIN_SIZE_SOFTCLIP;
                    break;

                case I:
                    insertTotal += element.getLength();
                    hasMinLengthIndel |= element.getLength() >= INDEL_TO_SC_MIN_SIZE_SOFTCLIP;
                    break;

                case M:
                    alignedTotal += element.getLength();
                    break;
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

    public static void buildIndelFrequencies(final Map<Integer,List<Read>> indelLengthReads, final Read read)
    {
        int maxIndelLength = read.indelCoords() != null ? read.indelCoords().Length : maxIndelLength(read.cigarElements());

        if(maxIndelLength >= INDEL_TO_SC_MIN_SIZE_SOFTCLIP)
        {
            List<Read> lengthReads = indelLengthReads.get(maxIndelLength);
            if(lengthReads == null)
            {
                lengthReads = Lists.newArrayList();
                indelLengthReads.put(maxIndelLength, lengthReads);
            }

            lengthReads.add(read);
        }
    }

    public static List<Read> findMaxFrequencyIndelReads(final Map<Integer,List<Read>> indelLengthReads)
    {
        if(indelLengthReads.isEmpty())
            return Collections.emptyList();

        int maxFrequency = 0;
        int indelLength = 0;

        for(Map.Entry<Integer,List<Read>> entry : indelLengthReads.entrySet())
        {
            if(entry.getValue().size() > maxFrequency)
            {
                indelLength = entry.getKey();
                maxFrequency = entry.getValue().size();
            }
        }

        return indelLengthReads.get(indelLength);
    }

    // OLD indel adjustment methods
    public static boolean convertEdgeIndelsToSoftClip(final Read read)
    {
        return convertEdgeIndelsToSoftClip(read, INDEL_TO_SC_MIN_SIZE_SOFTCLIP, INDEL_TO_SC_MAX_SIZE_SOFTCLIP);
    }

    public static boolean convertEdgeIndelsToSoftClip(final Read read, final int minIndelLength, final int maxIndelLength)
    {
        if(read.cigarElements().size() < 3)
            return false;

        int leftSoftClipLength = calcIndelToSoftClipLength(
                read.cigarElements().get(0), read.cigarElements().get(1), read.cigarElements().get(2),
                minIndelLength, maxIndelLength);

        int lastIndex = read.cigarElements().size() - 1;

        int rightSoftClipLength = calcIndelToSoftClipLength(
                read.cigarElements().get(lastIndex), read.cigarElements().get(lastIndex - 1), read.cigarElements().get(lastIndex - 2),
                minIndelLength, maxIndelLength);

        if(leftSoftClipLength > 0 || rightSoftClipLength > 0)
        {
            // read.setIndelUnclippedBounds(leftSoftClipLength, rightSoftClipLength);
            return true;
        }

        return false;
    }

    private static int calcIndelToSoftClipLength(
            final CigarElement edge, final CigarElement inside, final CigarElement next,
            final int minIndelLength, final int maxIndelLength)
    {
        if(edge.getOperator() != M)
            return 0;

        if(!inside.getOperator().isIndel())
            return 0;

        if(next.getOperator() != M)
            return 0;

        if(inside.getLength() < minIndelLength || inside.getLength() > maxIndelLength)
            return 0;

        return inside.getOperator() != D ? edge.getLength() + inside.getLength() : edge.getLength();
    }

    public static boolean hasDominantIndelReadAssembly(final JunctionAssembly assembly)
    {
        int indelReads = 0;
        int totalJuncReads = 0;

        for(SupportRead read : assembly.support())
        {
            if(read.type().isSplitSupport())
                ++totalJuncReads;

            if(read.type() == SupportType.INDEL)
                ++indelReads;
        }

        return indelReads >= 0.5 * totalJuncReads;
    }
}
