package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_MAP_QUALITY;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MAX_SUPPORT_FRAGMENT_DISTANCE;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DEFAULT_MAX_FRAGMENT_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_ALIGNMENT_BASES;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_INSERT_ALIGNMENT_OVERLAP;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_INSERT_LENGTH_SUPPORT;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_JUNCTION_SUPPORT;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_SOFT_CLIP_HIGH_QUAL_PERC;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_SUPPORTING_READ_DISTANCE;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.commons.cli.CommandLine;

public class ReadFilterConfig
{
    // initial read filtering to identify junction candidates
    public final int MinAlignmentBases;
    public final int MinMapQuality;
    public final int MinInsertAlignmentOverlap;
    public final int MinIndelLength;
    public final int MinSoftClipLength;
    public final double MinSoftClipHighQualPerc;

    // supporting read config
    public final int MinSupportingReadDistance;
    public final int MinInsertLengthSupport;

    // final junction filtering
    public final int MinJunctionSupport;

    private int mObservedFragLengthMin; // as set by the distribution min and max percentiles
    private int mObservedFragLengthMax;
    private int mMaxSupportingFragmentDistance;

    private static final String CFG_MIN_ALIGN_BASES = "min_align_bases";
    private static final String CFG_MIN_JUNCTION_FRAGS = "min_junction_frags";

    public ReadFilterConfig(
            final int minAlignmentBases, final int minMapQuality, final int minInsertAlignmentOverlap, final int minSoftClipLength,
            final double minSoftClipHighQualPerc, final int minSupportingReadDistance,
            final int minIndelLength, final int minJunctionSupport)
    {
        MinAlignmentBases = minAlignmentBases;
        MinMapQuality = minMapQuality;
        MinInsertAlignmentOverlap = minInsertAlignmentOverlap;
        MinSoftClipLength = minSoftClipLength;
        MinIndelLength = minIndelLength;
        MinSoftClipHighQualPerc = minSoftClipHighQualPerc;
        MinSupportingReadDistance = minSupportingReadDistance;

        MinInsertLengthSupport = MIN_INSERT_LENGTH_SUPPORT;
        MinJunctionSupport = minJunctionSupport;

        mMaxSupportingFragmentDistance = DEFAULT_MAX_FRAGMENT_LENGTH;
        mObservedFragLengthMax = DEFAULT_MAX_FRAGMENT_LENGTH;
        mObservedFragLengthMin = DEFAULT_READ_LENGTH;
    }

    public void setFragmentLengths(int minLength, int maxLength)
    {
        mObservedFragLengthMin = minLength;
        mObservedFragLengthMax = maxLength;
        mMaxSupportingFragmentDistance = min(mObservedFragLengthMax, MAX_SUPPORT_FRAGMENT_DISTANCE);
    }

    public int observedFragLengthMax() { return mObservedFragLengthMax; }
    public int observedFragLengthMin() { return mObservedFragLengthMin; }
    public int maxSupportingFragmentDistance() { return mMaxSupportingFragmentDistance; }

    public static ReadFilterConfig from(final ConfigBuilder configBuilder)
    {
        return new ReadFilterConfig(
                configBuilder.getInteger(CFG_MIN_ALIGN_BASES),
                MIN_MAP_QUALITY,
                MIN_INSERT_ALIGNMENT_OVERLAP,
                MIN_SOFT_CLIP_LENGTH,
                MIN_SOFT_CLIP_HIGH_QUAL_PERC,
                MIN_SUPPORTING_READ_DISTANCE,
                MIN_INDEL_LENGTH,
                configBuilder.getInteger(CFG_MIN_JUNCTION_FRAGS));
    }

    private static int configValue(final CommandLine cmd, final String config, final int defaultValue)
    {
        return cmd != null ? Integer.parseInt(cmd.getOptionValue(config, String.valueOf(defaultValue))) : defaultValue;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addInteger(CFG_MIN_ALIGN_BASES, "Junction fragment min aligned bases", MIN_ALIGNMENT_BASES);
        configBuilder.addInteger(CFG_MIN_JUNCTION_FRAGS, "Required fragments to keep a junction", MIN_JUNCTION_SUPPORT);
    }
}
