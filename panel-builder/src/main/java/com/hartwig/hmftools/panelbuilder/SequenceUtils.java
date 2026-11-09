package com.hartwig.hmftools.panelbuilder;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;
import static java.lang.System.arraycopy;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBasesInPlace;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionEndingAt;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionStartingAt;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class SequenceUtils
{
    public static SequenceData buildSequence(final RefGenomeInterface refGenome, final SequenceDefinition definition)
    {
        List<SequenceData> parts = definition.segments().stream()
                .map(segment -> buildSegmentSequence(refGenome, segment))
                .toList();

        byte[] sequence = new byte[definition.baseLength()];
        int sequenceIndex = 0;
        for(SequenceData part : parts)
        {
            arraycopy(part.bases(), 0, sequence, sequenceIndex, part.bases().length);
            sequenceIndex += part.bases().length;
        }

        boolean isNormal = parts.stream().allMatch(SequenceData::isNormal);
        int gcCount = parts.stream().mapToInt(SequenceData::gcCount).sum();
        return new SequenceData(sequence, isNormal, gcCount);
    }

    private static SequenceData buildSegmentSequence(final RefGenomeInterface refGenome, final SequenceSegment segment)
    {
        if(segment instanceof RefSegment refSegment)
        {
            SequenceData data = getSequence(refGenome, refSegment.region());
            if(refSegment.orientation() == Orientation.REVERSE)
            {
                // Safe to mutate in place: getSequence returns a fresh array. Reverse complement preserves GC count and normality.
                reverseComplementBasesInPlace(data.bases(), 0, data.bases().length);
            }
            return data;
        }
        else
        {
            return validateAndNormaliseSequence(((InsertSeqSegment) segment).sequence());
        }
    }

    private static SequenceData getSequence(final RefGenomeInterface refGenome, final ChrBaseRegion region)
    {
        byte[] sequence = refGenome.getBases(region.chromosome(), region.start(), region.end());
        if(sequence == null || sequence.length != region.baseLength())
        {
            throw new IllegalArgumentException("Attempt to create probe in unmapped region: " + region);
        }
        return validateAndNormaliseSequence(sequence);
    }

    public static SequenceData validateAndNormaliseSequence(final byte[] sequence)
    {
        // Combine these into one iteration to increase performance:
        //   - check if the sequence contains normal bases
        //   - convert to uppercase
        //   - count GC content

        byte[] normalised = new byte[sequence.length];
        boolean isNormal = true;
        int gcCount = 0;
        for(int i = 0; i < sequence.length; ++i)
        {
            byte b = sequence[i];
            byte n;
            switch(b)
            {
                case 'A':
                case 'T':
                    n = b;
                    break;
                case 'G':
                case 'C':
                    n = b;
                    ++gcCount;
                    break;
                case 'a':
                    n = 'A';
                    break;
                case 't':
                    n = 'T';
                    break;
                case 'g':
                    n = 'G';
                    ++gcCount;
                    break;
                case 'c':
                    n = 'C';
                    ++gcCount;
                    break;
                case 'N':
                    n = b;
                    isNormal = false;
                    break;
                case 'n':
                    n = 'N';
                    isNormal = false;
                    break;
                default:
                    n = (byte) Character.toUpperCase(b);
                    isNormal = false;
            }
            normalised[i] = n;
        }
        return new SequenceData(normalised, isNormal, gcCount);
    }

    public static SequenceData validateAndNormaliseSequence(final String sequence)
    {
        return validateAndNormaliseSequence(sequence.getBytes());
    }

    public static boolean isDnaSequenceNormal(final String sequence)
    {
        for(int i = 0; i < sequence.length(); ++i)
        {
            switch(sequence.charAt(i))
            {
                case 'A':
                case 'T':
                case 'G':
                case 'C':
                case 'a':
                case 't':
                case 'g':
                case 'c':
                    break;
                default:
                    return false;
            }
        }
        return true;
    }

    // Whether the probe sequence is within maxBasesDifference bases of a contiguous stretch of reference genome. True only for simple near-ref
    // probes (single region, single breakend + small insert, same-chromosome SNV/small INDEL/SV); anything more complex is dissimilar.
    public static boolean isSequenceSimilarToRef(final SequenceDefinition definition, int maxBasesDifference)
    {
        List<ChrBaseRegion> regions = definition.regions();
        int regionBases = regions.stream().mapToInt(ChrBaseRegion::baseLength).sum();
        int insertBases = definition.baseLength() - regionBases;

        if(regions.size() == 1)
        {
            // Exact single region (no insert) or single breakend (region + insert): difference is just the insert.
            return insertBases <= maxBasesDifference;
        }
        if(regions.size() == 2)
        {
            ChrBaseRegion a = regions.get(0);
            ChrBaseRegion b = regions.get(1);
            if(!a.chromosome().equals(b.chromosome()))
            {
                // Translocation.
                return false;
            }
            // Clamp to >=0 because theoretically the regions could overlap in the case of an SV.
            ChrBaseRegion lower = a.start() <= b.start() ? a : b;
            ChrBaseRegion upper = a.start() <= b.start() ? b : a;
            int deleteBases = max(upper.start() - lower.end() - 1, 0);
            return max(insertBases, deleteBases) <= maxBasesDifference;
        }
        // More than two regions (spliced) or otherwise complex: dissimilar.
        return false;
    }

    public static SequenceDefinition buildIndelProbe(final String chromosome, int position, final String ref, final String alt,
            int probeLength)
    {
        int altLength = alt.length();
        int refLength = ref.length();
        int startLength = probeLength / 2 - altLength / 2;

        ChrBaseRegion startRegion = regionEndingAt(chromosome, position - 1, startLength);
        int endBaseLength = probeLength - startRegion.baseLength() - altLength;

        int postPosition = position + refLength;
        ChrBaseRegion endRegion = regionStartingAt(chromosome, postPosition, endBaseLength);

        SequenceDefinition definition = SequenceDefinition.variant(startRegion, Orientation.FORWARD, alt, endRegion, Orientation.FORWARD);

        if(definition.baseLength() != probeLength)
        {
            throw new IllegalArgumentException(format(
                    "variant(%s:%d %s->%s) invalid sequenceLength(%d)",
                    chromosome, position, ref, alt, definition.baseLength()));
        }

        return definition;
    }

    public static SequenceDefinition buildSvProbe(final String chrStart, int positionStart, Orientation orientStart, final String chrEnd,
            int positionEnd, Orientation orientEnd, final String insertSequence, int probeLength)
    {
        int halfProbeLength = probeLength / 2;
        int insSeqLength = insertSequence.length();
        int halfInsSeqLength = insSeqLength / 2;
        int halfNonInsSeqLength = halfProbeLength - halfInsSeqLength;

        // +1 to -1 - do nothing
        // -1 to +1 - a DUP

        // +1 to +1 - start normal, add insert and then reverse compliment of the other side
        // -1 to -1 - alt insert sequence then ref

        ChrBaseRegion startRegion;
        Orientation startOrient;
        ChrBaseRegion endRegion;
        Orientation endOrient;

        if(orientStart == Orientation.FORWARD)
        {
            startRegion = regionEndingAt(chrStart, positionStart, halfNonInsSeqLength);
            startOrient = Orientation.FORWARD;

            int endBaseLength = probeLength - startRegion.baseLength() - insSeqLength;

            if(orientEnd == Orientation.REVERSE)
            {
                endRegion = regionStartingAt(chrEnd, positionEnd, endBaseLength);
                endOrient = Orientation.FORWARD;
            }
            else
            {
                endRegion = regionEndingAt(chrEnd, positionEnd, endBaseLength);
                endOrient = Orientation.REVERSE;
            }
        }
        else
        {
            if(orientEnd == Orientation.FORWARD)
            {
                // swap ends and treat as +1/-1
                startRegion = regionEndingAt(chrEnd, positionEnd, halfNonInsSeqLength);
                startOrient = Orientation.FORWARD;
            }
            else
            {
                // -1/-1 - start with the reversed bases from the end breakend
                startRegion = regionStartingAt(chrEnd, positionEnd, halfNonInsSeqLength);
                startOrient = Orientation.REVERSE;
            }
            int endBaseLength = probeLength - startRegion.baseLength() - insSeqLength;
            endRegion = regionStartingAt(chrStart, positionStart, endBaseLength);
            endOrient = Orientation.FORWARD;
        }

        SequenceDefinition definition = SequenceDefinition.variant(startRegion, startOrient, insertSequence, endRegion, endOrient);

        if(definition.baseLength() != probeLength)
        {
            throw new IllegalArgumentException("Invalid probe");
        }

        return definition;
    }

    public static SequenceDefinition buildSglProbe(final String chromosome, int position, Orientation orientation,
            final String insertSequence, int probeLength)
    {
        int halfProbeLength = probeLength / 2;
        int insSeqLength = min(insertSequence.length(), halfProbeLength);
        int refBaseLength = probeLength - insSeqLength;

        SequenceDefinition definition;

        if(orientation == Orientation.FORWARD)
        {
            ChrBaseRegion startRegion = regionEndingAt(chromosome, position, refBaseLength);
            String insert = insertSequence.substring(0, insSeqLength);
            if(insert.isEmpty())
            {
                definition = SequenceDefinition.singleRegion(startRegion);
            }
            else
            {
                definition = SequenceDefinition.forwardSgl(startRegion, insert);
            }
        }
        else
        {
            String insert = insertSequence.substring(insertSequence.length() - insSeqLength);
            ChrBaseRegion endRegion = regionStartingAt(chromosome, position, refBaseLength);
            if(insert.isEmpty())
            {
                definition = SequenceDefinition.singleRegion(endRegion);
            }
            else
            {
                definition = SequenceDefinition.reverseSgl(insert, endRegion);
            }
        }

        if(definition.baseLength() != probeLength)
        {
            throw new IllegalArgumentException("Invalid probe");
        }

        return definition;
    }
}
