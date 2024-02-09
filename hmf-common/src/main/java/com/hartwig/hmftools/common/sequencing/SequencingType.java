package com.hartwig.hmftools.common.sequencing;

import java.util.Arrays;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public enum SequencingType
{
    ILLUMINA,
    ULTIMA;

    public static final String SEQUENCING_TYPE_CFG = "sequencing_type";

    public static final String SEQUENCING_TYPE_DESC_CFG = "Sequencing type: {}"
            + Arrays.stream(SequencingType.values()).map(x -> x.toString()).collect(Collectors.joining(", "));

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SEQUENCING_TYPE_CFG, false, SEQUENCING_TYPE_DESC_CFG, ILLUMINA.toString());
    }

}
