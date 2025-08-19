package com.hartwig.hmftools.lilac.variant;

import com.hartwig.hmftools.lilac.hla.HlaGene;

public class CopyNumberData
{
    public final HlaGene Gene;
    public final double MinCopyNumber;
    public final double MinMinorAlleleCopyNumber;

    public CopyNumberData(final HlaGene gene, final double minCopyNumber, final double minMinorAlleleCopyNumber)
    {
        Gene = gene;
        MinCopyNumber = minCopyNumber;
        MinMinorAlleleCopyNumber = minMinorAlleleCopyNumber;
    }

    @Override
    public String toString()
    {
        return String.format("gene(%s) min(%.3f) minor(%.3f)", Gene, MinCopyNumber, MinMinorAlleleCopyNumber);
    }
}
