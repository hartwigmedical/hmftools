package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MAX_SUPPORT_FRAGMENT_DISTANCE;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DEFAULT_MAX_FRAGMENT_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_ALIGNMENT_BASES;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_INSERT_ALIGNMENT_OVERLAP;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_INSERT_LENGTH_SUPPORT;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_JUNCTION_SUPPORT;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_MAP_QUALITY;
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
    public final int MinSoftClipHighQual;
    public final double MinSoftClipHighQualPerc;

    // supporting read config
    public final int MinSupportingReadDistance;
    public final int MinInsertLengthSupport;

    // final junction filtering
    public final int MinJunctionSupport;

    private int mFragmentLengthMin; // as set by the distribution min and max percentiles
    private int mFragmentLengthMax;
    private int mMaxSupportingFragmentDistance;

    private static final String CFG_MIN_ALIGN_BASES = "min_align_bases";
    private static final String CFG_MIN_JUNCTION_FRAGS = "min_junction_frags";

    public ReadFilterConfig(
            final int minAlignmentBases, final int minMapQuality, final int minInsertAlignmentOverlap, final int minSoftClipLength,
            final int minSoftClipHighQual, final double minSoftClipHighQualPerc, final int minSupportingReadDistance,
            final int minIndelLength, final int minJunctionSupport)
    {
        MinAlignmentBases = minAlignmentBases;
        MinMapQuality = minMapQuality;
        MinInsertAlignmentOverlap = minInsertAlignmentOverlap;
        MinSoftClipLength = minSoftClipLength;
        MinIndelLength = minIndelLength;
        MinSoftClipHighQual = minSoftClipHighQual;
        MinSoftClipHighQualPerc = minSoftClipHighQualPerc;
        MinSupportingReadDistance = minSupportingReadDistance;

        MinInsertLengthSupport = MIN_INSERT_LENGTH_SUPPORT;
        MinJunctionSupport = minJunctionSupport;

        mMaxSupportingFragmentDistance = DEFAULT_MAX_FRAGMENT_LENGTH;
        mFragmentLengthMax = DEFAULT_MAX_FRAGMENT_LENGTH;
        mFragmentLengthMin = DEFAULT_READ_LENGTH;
    }

    public void setFragmentLengths(int minLength, int maxLength)
    {
        mFragmentLengthMin = minLength;
        mFragmentLengthMax = maxLength;
        mMaxSupportingFragmentDistance = min(mFragmentLengthMax, MAX_SUPPORT_FRAGMENT_DISTANCE);
    }

    public int fragmentLengthMax() { return mFragmentLengthMax; }
    public int fragmentLengthMin() { return mFragmentLengthMin; }
    public int maxSupportingFragmentDistance() { return mMaxSupportingFragmentDistance; }

    public static ReadFilterConfig from(final ConfigBuilder configBuilder)
    {
        return new ReadFilterConfig(
                configBuilder.getInteger(CFG_MIN_ALIGN_BASES),
                MIN_MAP_QUALITY,
                MIN_INSERT_ALIGNMENT_OVERLAP,
                MIN_SOFT_CLIP_LENGTH,
                LOW_BASE_QUAL_THRESHOLD,
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
