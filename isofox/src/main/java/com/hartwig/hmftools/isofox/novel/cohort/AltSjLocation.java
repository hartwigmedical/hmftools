package com.hartwig.hmftools.isofox.novel.cohort;

import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunctionType;

public class AltSjLocation
{
    public final String GeneId;
    public final BaseRegion Location;
    public final AltSpliceJunctionType Type;
    public final String Key;

    public AltSjLocation(final BaseRegion location, final String geneId, final AltSpliceJunctionType type)
    {
        GeneId = geneId;
        Location = location;
        Type = type;

        Key = String.format("%s-%d-%d", location.Chromosome, location.start(), location.end());
    }

    public static AltSjLocation fromCsv(final String[] data, int geneIndex, int chrIndex, int posStartIndex, int posEndIndex, int typeIndex)
    {
        return new AltSjLocation(
                new BaseRegion(data[chrIndex], Integer.parseInt(data[posStartIndex]), Integer.parseInt(data[posEndIndex])),
                data[geneIndex], AltSpliceJunctionType.valueOf(data[typeIndex]));
    }
}
