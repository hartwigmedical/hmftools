package com.hartwig.hmftools.bamtools.synthetic;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.Collection;
import java.util.stream.Collectors;

public enum RegionType
{
    PANEL,
    DRIVER_GENE,
    COPY_NUMBER,
    CN_BACKBONE,
    GERMLINE_HET_SITE,
    SOMATIC_MUTATION,
    GERMLINE_MUTATION,
    SOMATIC_SV,
    FUSION,
    DISRUPTION,
    GERMLINE_SV;

    public static String typesStr(final Collection<RegionType> types)
    {
        return types.stream().map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM));
    }
}
