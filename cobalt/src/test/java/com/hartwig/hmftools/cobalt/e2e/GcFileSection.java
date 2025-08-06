package com.hartwig.hmftools.cobalt.e2e;

public interface GcFileSection extends ChromosomeSection
{
    @Override
    default String line(final int position)
    {
        return String.format("chr%s\t%d\t%.2f\t1\t1", chromosome().toString(), position, gcRatio(position));
    }

    double gcRatio(int position);
}
