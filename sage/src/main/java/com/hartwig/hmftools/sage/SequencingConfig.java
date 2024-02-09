package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.common.sequencing.SequencingType.SEQUENCING_TYPE_CFG;

import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class SequencingConfig
{
    public final boolean HasUMIs;
    public final SequencingType Type;

    private static final String TRACK_UMIS = "track_umis";

    public SequencingConfig(final boolean hasUMIs, final SequencingType type)
    {
        HasUMIs = hasUMIs;
        Type = type;
    }

    public static SequencingConfig from(final ConfigBuilder configBuilder)
    {
        SequencingType sequencingType = SequencingType.valueOf(configBuilder.getValue(SEQUENCING_TYPE_CFG));
        boolean useUMIs = configBuilder.hasFlag(TRACK_UMIS);

        return new SequencingConfig(useUMIs, sequencingType);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addFlag(TRACK_UMIS, "Record counts of UMI types");
        SequencingType.registerConfig(configBuilder);
    }
}
