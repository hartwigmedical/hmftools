package com.hartwig.hmftools.lilac.variant;

import com.sun.tools.javah.Gen;

public class CopyNumberData
{
    public final String Gene;
    public final double MinCopyNumber;
    public final double MinMinorAlleleCopyNumber;

    public CopyNumberData(final String gene, final double minCopyNumber, final double minMinorAlleleCopyNumber)
    {
        Gene = gene;
        MinCopyNumber = minCopyNumber;
        MinMinorAlleleCopyNumber = minMinorAlleleCopyNumber;
    }

    public String toString()
    {
        return String.format("gene(%s) min(%.3f) minor(%.3f)", Gene, MinCopyNumber, MinMinorAlleleCopyNumber);
    }
}
