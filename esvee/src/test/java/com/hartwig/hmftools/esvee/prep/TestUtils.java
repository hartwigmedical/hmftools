package com.hartwig.hmftools.esvee.prep;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.prep.types.ReadFilterConfig;

public final class TestUtils
{
    public static final int PARTITION_SIZE = 10000;
    public static final ChrBaseRegion REGION_1 = new ChrBaseRegion(CHR_1, 1, PARTITION_SIZE - 1);
    public static final ChrBaseRegion REGION_2 = new ChrBaseRegion(CHR_1, REGION_1.end() + 1, REGION_1.end() + PARTITION_SIZE);
    public static final ChrBaseRegion REGION_3 = new ChrBaseRegion(CHR_1, REGION_2.end() + 1, REGION_2.end() + PARTITION_SIZE);

    public static final ConfigBuilder READ_FILTERS_CONFIG = new ConfigBuilder();

    static
    {
        ReadFilterConfig.addConfig(READ_FILTERS_CONFIG);
    }

    public static final ReadFilterConfig READ_FILTERS = ReadFilterConfig.from(READ_FILTERS_CONFIG);
    public static final HotspotCache HOTSPOT_CACHE = new HotspotCache(null);
    public static final BlacklistLocations BLACKLIST_LOCATIONS = new BlacklistLocations(null);
}
