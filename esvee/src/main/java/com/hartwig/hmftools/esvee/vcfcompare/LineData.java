package com.hartwig.hmftools.esvee.vcfcompare;

import com.hartwig.hmftools.esvee.vcfcompare.line.LineLink;

public class LineData
{
    public final boolean HasPolyAT;
    public LineLink LinkedLineBreakends;
    public LineLink InferredLinkedLineBreakends;

    public LineData(final boolean hasPolyAT)
    {
        HasPolyAT = hasPolyAT;
        LinkedLineBreakends = null;
        InferredLinkedLineBreakends = null;
    }

    public boolean hasLineLink() { return LinkedLineBreakends != null; }
    public boolean hasInferredLineLink() { return InferredLinkedLineBreakends != null; }
}
