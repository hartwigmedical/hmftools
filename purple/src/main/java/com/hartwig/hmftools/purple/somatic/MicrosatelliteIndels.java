package com.hartwig.hmftools.purple.somatic;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.variant.SageVcfTags.REPEAT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.REPEAT_SEQUENCE;
import static com.hartwig.hmftools.purple.PurpleConstants.TUMOR_MSI_LOAD_MIN_VAF;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeConstants.MB_PER_GENOME;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.redux.MsiModelPrediction;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.purple.targeted.TargetRegionsData;

import org.apache.logging.log4j.util.Strings;

public class MicrosatelliteIndels
{
    private final TargetRegionsData mTargetRegions;
    private final MsiModelPrediction mMsiModelPrediction;
    private int mIndelCount;

    private static final int MIN_SEQUENCE_LENGTH_FOR_LONG_REPEATS = 2;
    private static final int MAX_SEQUENCE_LENGTH_FOR_LONG_REPEATS = 4;
    private static final int MIN_REPEAT_COUNT_FOR_LONG_REPEATS = 4;
    private static final int MIN_REPEAT_COUNT_FOR_SHORT_REPEATS = 5;

    private static final int MAX_REF_ALT_LENGTH = 50;

    public MicrosatelliteIndels(final TargetRegionsData referenceData, final MsiModelPrediction msiModelPrediction)
    {
        mTargetRegions = referenceData;
        mMsiModelPrediction = msiModelPrediction;
        mIndelCount = 0;
    }

    public int msiIndelCount() { return mIndelCount; }

    public double calcMsiIndelsPerMb()
    {
        if(mTargetRegions.hasTargetRegions())
        {
            return mMsiModelPrediction != null ? mMsiModelPrediction.PredictedMsiIndelsPerMb : 0;
        }
        else
        {
            return mIndelCount / MB_PER_GENOME;
        }
    }

    public void processVariant(final SomaticVariant variant)
    {
        if(variant.type() != VariantType.INDEL)
            return;

        int altLength = variant.decorator().alt().length();
        int refLength = variant.decorator().ref().length();

        if(refLength >= MAX_REF_ALT_LENGTH || altLength >= MAX_REF_ALT_LENGTH)
            return;

        double vaf = variant.alleleFrequency();

        if(vaf < TUMOR_MSI_LOAD_MIN_VAF)
            return;

        int repeatCount = variant.context().getAttributeAsInt(REPEAT_COUNT, 0);
        int repeatSequenceLength = variant.context().getAttributeAsString(REPEAT_SEQUENCE, Strings.EMPTY).length();

        if(!repeatContextIsRelevant(repeatCount, repeatSequenceLength))
            return;

        if(mTargetRegions.hasTargetRegions() && PPL_LOGGER.isTraceEnabled())
        {
            PPL_LOGGER.trace(format("indel(%s) af(%.2f) included in target-regions TMB", variant, vaf));
        }

        mIndelCount++;
    }

    private static boolean repeatContextIsRelevant(int repeatCount, int repeatSequenceLength)
    {
        boolean longRepeatRelevant = repeatSequenceLength >= MIN_SEQUENCE_LENGTH_FOR_LONG_REPEATS
                && repeatSequenceLength <= MAX_SEQUENCE_LENGTH_FOR_LONG_REPEATS
                && repeatCount >= MIN_REPEAT_COUNT_FOR_LONG_REPEATS;

        boolean shortRepeatRelevant = repeatSequenceLength == 1 && repeatCount >= MIN_REPEAT_COUNT_FOR_SHORT_REPEATS;

        return longRepeatRelevant | shortRepeatRelevant;
    }

    @VisibleForTesting
    static boolean repeatContextIsRelevant(int repeatCount, String sequence)
    {
        return repeatContextIsRelevant(repeatCount, sequence.length());
    }
}
