package com.hartwig.hmftools.geneutils.panelfinder;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class PanelData extends ChrBaseRegion
{
    public final String Label;

    public PanelData(final ChrBaseRegion region, final String label)
    {
        super(region.Chromosome, region.start(), region.end());

        Label = label;
    }

    public String toString() { return format("%s %s", super.toString(), Label); }

    public static String toString(final List<PanelData> regions)
    {
        return regions.stream().map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM));
    }

}
