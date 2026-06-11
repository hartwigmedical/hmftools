package com.hartwig.hmftools.esvee.common.saga;

import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_JUNCTION_INDEL_DISTANCE;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_JUNCTION_INDEL_MAX_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_JUNCTION_OVERLAP_MIN;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_JUNCTION_OVERLAP_MIN_LOWER;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_LENGTH_MIN_RATIO;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_SCORE_MIN_BASELINE;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_SCORE_MIN_RATIO;

public record SagaSequenceMatcherConfig(
        double alignLengthRatioMin,
        int alignScoreMin,
        double alignScoreRatioMin,
        int junctionOverlapMin,
        int junctionOverlapMinLower,
        int junctionIndelDistance,
        int junctionIndelLengthMax
)
{
    public SagaSequenceMatcherConfig
    {
        if(!(alignLengthRatioMin >= 0 && alignLengthRatioMin <= 1))
        {
            throw new IllegalArgumentException();
        }
        if(alignScoreMin < 0)
        {
            throw new IllegalArgumentException();
        }
        if(!(alignScoreRatioMin >= 0 && alignScoreRatioMin <= 1))
        {
            throw new IllegalArgumentException();
        }
        if(junctionIndelDistance < 0)
        {
            throw new IllegalArgumentException();
        }
        if(junctionIndelLengthMax < 0)
        {
            throw new IllegalArgumentException();
        }
    }

    static SagaSequenceMatcherConfig DEFAULT = new SagaSequenceMatcherConfig(
            SAGA_ALIGN_LENGTH_MIN_RATIO,
            SAGA_ALIGN_SCORE_MIN_BASELINE,
            SAGA_ALIGN_SCORE_MIN_RATIO,
            SAGA_ALIGN_JUNCTION_OVERLAP_MIN,
            SAGA_ALIGN_JUNCTION_OVERLAP_MIN_LOWER,
            SAGA_ALIGN_JUNCTION_INDEL_DISTANCE,
            SAGA_ALIGN_JUNCTION_INDEL_MAX_LENGTH
    );
}
