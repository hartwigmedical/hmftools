package com.hartwig.hmftools.svprep.reads;

import static java.lang.Math.min;

import static com.hartwig.hmftools.svprep.SvConstants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.svprep.SvConstants.MAX_SUPPORT_FRAGMENT_DISTANCE;
import static com.hartwig.hmftools.svprep.SvConstants.DEFAULT_MAX_FRAGMENT_LENGTH;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_ALIGNMENT_BASES;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_INDEL_LENGTH;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_INSERT_ALIGNMENT_OVERLAP;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_INSERT_LENGTH_SUPPORT;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_JUNCTION_SUPPORT;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_MAP_QUALITY;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_SOFT_CLIP_HIGH_QUAL_PERC;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_SOFT_CLIP_MIN_BASE_QUAL;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_SUPPORTING_READ_DISTANCE;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

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

    public static ReadFilterConfig from(final CommandLine cmd)
    {
        return new ReadFilterConfig(
                configValue(cmd, CFG_MIN_ALIGN_BASES, MIN_ALIGNMENT_BASES),
                MIN_MAP_QUALITY,
                MIN_INSERT_ALIGNMENT_OVERLAP,
                MIN_SOFT_CLIP_LENGTH,
                MIN_SOFT_CLIP_MIN_BASE_QUAL,
                MIN_SOFT_CLIP_HIGH_QUAL_PERC,
                MIN_SUPPORTING_READ_DISTANCE,
                MIN_INDEL_LENGTH,
                configValue(cmd, CFG_MIN_JUNCTION_FRAGS, MIN_JUNCTION_SUPPORT));
    }

    private static int configValue(final CommandLine cmd, final String config, final int defaultValue)
    {
        return cmd != null ? Integer.parseInt(cmd.getOptionValue(config, String.valueOf(defaultValue))) : defaultValue;
    }

    public static void addCmdLineArgs(final Options options)
    {
        options.addOption(
                CFG_MIN_ALIGN_BASES, true, "Junction fragment min aligned bases, default: " + MIN_ALIGNMENT_BASES);

        options.addOption(
                CFG_MIN_JUNCTION_FRAGS, true, "Required fragments to keep a junction, default: " + MIN_JUNCTION_SUPPORT);
    }
}
