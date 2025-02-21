package com.hartwig.hmftools.common.basequal.jitter;

import static com.hartwig.hmftools.common.sequencing.SequencingType.BIOMODAL;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ILLUMINA;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;

import java.util.Collections;
import java.util.EnumSet;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.sequencing.SequencingType;

public enum ConsensusType
{
    IGNORE,
    NONE,
    SINGLE,
    DUAL,
    HIGH_QUAL;

    public static EnumSet<ConsensusType> consensusTypes(final JitterAnalyserConfig config)
    {
        SequencingType sequencingType = config.Sequencing;
        EnumSet<ConsensusType> consensusTypes = Sets.newEnumSet(Collections.emptyList(), ConsensusType.class);
        consensusTypes.add(NONE);
        if(sequencingType == ILLUMINA && config.UsesDuplexUMIs)
        {
            consensusTypes.add(SINGLE);
            consensusTypes.add(DUAL);
        }
        else if(sequencingType == SBX)
        {
            consensusTypes.add(DUAL);
        }
        else if(sequencingType == BIOMODAL)
        {
            consensusTypes.add(HIGH_QUAL);
        }

        return consensusTypes;
    }
}
