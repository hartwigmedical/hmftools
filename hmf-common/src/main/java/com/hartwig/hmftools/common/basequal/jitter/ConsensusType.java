package com.hartwig.hmftools.common.basequal.jitter;

import static com.hartwig.hmftools.common.sequencing.SequencingType.BIOMODAL;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;

import java.util.Arrays;
import java.util.Collections;
import java.util.EnumSet;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.sequencing.SequencingType;

public enum ConsensusType
{
    IGNORE(),
    NONE(SequencingType.values()),
    DUPLEX(SBX),
    HIGH_QUAL(BIOMODAL);

    public final EnumSet<SequencingType> SequencingTypes;

    ConsensusType(final SequencingType... sequencingTypes)
    {
        SequencingTypes = Sets.newEnumSet(Arrays.asList(sequencingTypes), SequencingType.class);
    }

    public static EnumSet<ConsensusType> consensusTypesFromSequencing(final SequencingType sequencingType)
    {
        EnumSet<ConsensusType> consensusTypes = Sets.newEnumSet(Collections.emptyList(), ConsensusType.class);
        for(ConsensusType consensusType : ConsensusType.values())
        {
            if(consensusType.SequencingTypes.contains(sequencingType))
                consensusTypes.add(consensusType);
        }

        return consensusTypes;
    }
}
