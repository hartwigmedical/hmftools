package com.hartwig.hmftools.common.sequencing;

import static com.hartwig.hmftools.common.utils.config.ConfigItem.enumValuesAsStr;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public enum SequencingType
{
    ILLUMINA,
    ULTIMA,
    SBX;

    public static final String SEQUENCING_TYPE_CFG = "sequencing_type";

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(
                SEQUENCING_TYPE_CFG, false, enumValuesAsStr(SequencingType.values(), "Sequencing types"),
                ILLUMINA.toString());
    }

    public static SequencingType parseConfig(final ConfigBuilder configBuilder)
    {
        if(configBuilder.hasValue(SEQUENCING_TYPE_CFG))
            return SequencingType.valueOf(configBuilder.getValue(SEQUENCING_TYPE_CFG));

        return ILLUMINA;
    }
}
