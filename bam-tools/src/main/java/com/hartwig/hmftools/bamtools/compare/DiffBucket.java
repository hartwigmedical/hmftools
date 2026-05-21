package com.hartwig.hmftools.bamtools.compare;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public enum DiffBucket
{
    ORIG_ONLY_PRIMARY,
    ORIG_ONLY_SUPP,
    ORIG_ONLY_SEC,
    NEW_ONLY_PRIMARY,
    NEW_ONLY_SUPP,
    NEW_ONLY_SEC,

    RESCUE_ORIG_UNMAPPED,
    RESCUE_NEW_UNMAPPED,
    CIGAR_MATCH_DIFF_LOCUS,
    SPLICE_VS_CONTIG,
    JUNCTION_SHIFT,
    PLACEMENT,
    SA_DIFF,
    MATE_ECHO,
    ATTR_ONLY,
    MIXED;

    public static DiffBucket classify(
            final MismatchType mismatchType, final SAMRecord origRead, final SAMRecord newRead, final List<String> diffs)
    {
        if(mismatchType == MismatchType.ORIG_ONLY)
        {
            if(origRead.isSecondaryAlignment())
                return ORIG_ONLY_SEC;
            if(origRead.getSupplementaryAlignmentFlag())
                return ORIG_ONLY_SUPP;
            return ORIG_ONLY_PRIMARY;
        }

        if(mismatchType == MismatchType.NEW_ONLY)
        {
            if(origRead.isSecondaryAlignment())
                return NEW_ONLY_SEC;
            if(origRead.getSupplementaryAlignmentFlag())
                return NEW_ONLY_SUPP;
            return NEW_ONLY_PRIMARY;
        }

        Set<String> fields = diffFieldNames(diffs);

        if(fields.contains("unmapped"))
        {
            if(origRead.getReadUnmappedFlag())
                return RESCUE_ORIG_UNMAPPED;
            if(newRead.getReadUnmappedFlag())
                return RESCUE_NEW_UNMAPPED;
        }

        // bucket on the non-mapQuality fields so mapQuality + a structural diff doesn't fall through to MIXED.
        Set<String> structuralFields = new HashSet<>(fields);
        structuralFields.remove("mapQuality");

        boolean bothMapped = !origRead.getReadUnmappedFlag() && !newRead.getReadUnmappedFlag();

        if(bothMapped && structuralFields.contains("cigar"))
        {
            int origN = countNBlocks(origRead.getCigar());
            int newN = countNBlocks(newRead.getCigar());
            if(origN != newN)
                return SPLICE_VS_CONTIG;
            if(origN > 0 && !junctionPositionsMatch(origRead, newRead))
                return JUNCTION_SHIFT;
        }

        if(structuralFields.contains("coords"))
        {
            // same alignment shape at a different locus = paralog pick / cross-locus disagreement.
            // Require both sides mapped and matching cigar strings; not an unmapped-flag diff in disguise.
            if(!structuralFields.contains("cigar") && !fields.contains("unmapped")
                    && bothMapped
                    && origRead.getCigarString().equals(newRead.getCigarString()))
                return CIGAR_MATCH_DIFF_LOCUS;
            return PLACEMENT;
        }

        if(structuralFields.isEmpty())
            return MIXED;

        if(structuralFields.stream().allMatch(DiffBucket::isMateField))
            return MATE_ECHO;

        if(structuralFields.size() == 1 && structuralFields.contains("attrib_SA"))
            return SA_DIFF;

        if(structuralFields.stream().allMatch(f -> f.startsWith("attrib_")))
            return ATTR_ONLY;

        return MIXED;
    }

    private static Set<String> diffFieldNames(final List<String> diffs)
    {
        Set<String> fields = new HashSet<>();
        for(String entry : diffs)
        {
            int paren = entry.indexOf('(');
            if(paren > 0)
                fields.add(entry.substring(0, paren));
        }
        return fields;
    }

    private static boolean isMateField(final String field)
    {
        return field.equals("matePos") || field.equals("mateChr") || field.equals("mateNegStrand")
                || field.equals("attrib_MC") || field.equals("mateUnmapped");
    }

    private static int countNBlocks(final Cigar cigar)
    {
        int count = 0;
        for(CigarElement e : cigar.getCigarElements())
        {
            if(e.getOperator() == CigarOperator.N)
                ++count;
        }
        return count;
    }

    private static boolean junctionPositionsMatch(final SAMRecord orig, final SAMRecord newRead)
    {
        return encodeJunctionIntervals(orig).equals(encodeJunctionIntervals(newRead));
    }

    private static String encodeJunctionIntervals(final SAMRecord read)
    {
        StringBuilder sb = new StringBuilder();
        int pos = read.getAlignmentStart();
        for(CigarElement e : read.getCigar().getCigarElements())
        {
            if(e.getOperator() == CigarOperator.N)
            {
                sb.append(pos).append('-').append(pos + e.getLength() - 1).append(';');
                pos += e.getLength();
            }
            else if(e.getOperator().consumesReferenceBases())
            {
                pos += e.getLength();
            }
        }
        return sb.toString();
    }
}
