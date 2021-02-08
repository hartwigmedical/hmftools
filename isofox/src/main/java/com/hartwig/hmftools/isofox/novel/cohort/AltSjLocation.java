package com.hartwig.hmftools.isofox.novel.cohort;

import com.hartwig.hmftools.common.utils.sv.BaseRegion;

public class AltSjLocation
{
    public final String GeneId;
    public final BaseRegion Location;
    public final String Key;

    public AltSjLocation(final BaseRegion location, final String geneId)
    {
        GeneId = geneId;
        Location = location;

        Key = String.format("%s-%d-%d", location.Chromosome, location.start(), location.end());
    }

    public static AltSjLocation fromCsv(final String[] data, int geneIndex, int chrIndex, int posStartIndex, int posEndIndex)
    {
        return new AltSjLocation(
                new BaseRegion(data[chrIndex], Integer.parseInt(data[posStartIndex]), Integer.parseInt(data[posEndIndex])),
                data[geneIndex]);
    }
}
