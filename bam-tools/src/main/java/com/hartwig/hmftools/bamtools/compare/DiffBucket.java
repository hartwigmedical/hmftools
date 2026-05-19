package com.hartwig.hmftools.bamtools.compare;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

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
    PLACEMENT,
    MATE_ECHO,
    FLAG_ONLY,
    ATTR_ONLY,
    MIXED;

    public static DiffBucket classify(
            final MismatchType type, final SAMRecord orig, final SAMRecord newRead, final List<String> diffs)
    {
        if(type == MismatchType.ORIG_ONLY)
        {
            if(orig.isSecondaryAlignment())
                return ORIG_ONLY_SEC;
            if(orig.getSupplementaryAlignmentFlag())
                return ORIG_ONLY_SUPP;
            return ORIG_ONLY_PRIMARY;
        }

        if(type == MismatchType.NEW_ONLY)
        {
            if(orig.isSecondaryAlignment())
                return NEW_ONLY_SEC;
            if(orig.getSupplementaryAlignmentFlag())
                return NEW_ONLY_SUPP;
            return NEW_ONLY_PRIMARY;
        }

        Set<String> fields = diffFieldNames(diffs);

        if(fields.contains("unmapped"))
        {
            if(orig.getReadUnmappedFlag())
                return RESCUE_ORIG_UNMAPPED;
            if(newRead.getReadUnmappedFlag())
                return RESCUE_NEW_UNMAPPED;
        }

        if(fields.contains("coords"))
            return PLACEMENT;

        if(fields.isEmpty())
            return MIXED;

        if(fields.stream().allMatch(DiffBucket::isMateField))
            return MATE_ECHO;

        if(fields.stream().allMatch(DiffBucket::isFlagField))
            return FLAG_ONLY;

        if(fields.stream().allMatch(f -> f.startsWith("attrib_")))
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

    private static boolean isFlagField(final String field)
    {
        return field.equals("negStrand") || field.equals("duplicate") || field.equals("unmapped");
    }
}
