package com.hartwig.hmftools.cobalt.e2e;

public interface GcFileSection extends ChromosomeSection
{
    @Override
    default String line(final int position)
    {
        return String.format("chr%s\t%d\t%.2f\t1\t%.2f", chromosome().toString(), position, gcRatio(position), mappablePercentage(position));
    }

    double gcRatio(int position);

    default double mappablePercentage(int position)
    {
        return 1.0;
    }
}
