package com.hartwig.hmftools.svprep.reads;

import static com.hartwig.hmftools.svprep.SvConstants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.svprep.SvConstants.JUNCTION_SUPPORT_CAP;
import static com.hartwig.hmftools.svprep.SvConstants.MAX_DISCORDANT_READ_DISTANCE;
import static com.hartwig.hmftools.svprep.SvConstants.MAX_FRAGMENT_LENGTH;
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
    public final int MaxDiscordantFragmentDistance;

    // final junction filtering
    public final int MinJunctionSupport;
    public final int MaxJunctionSupportingReads;

    private int mFragmentLengthMin;
    private int mFragmentLengthMax;

    private static final String CFG_MIN_ALIGN_BASES = "min_align_bases";
    private static final String CFG_MAX_JUNCTION_SUPPORTING_READS = "max_junction_support_reads";
    private static final String CFG_MIN_JUNCTION_FRAGS = "min_junction_frags";

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
                configValue(cmd, CFG_MIN_JUNCTION_FRAGS, MIN_JUNCTION_SUPPORT),
                configValue(cmd, CFG_MAX_JUNCTION_SUPPORTING_READS, JUNCTION_SUPPORT_CAP),
                MAX_DISCORDANT_READ_DISTANCE);
    }

    private static int configValue(final CommandLine cmd, final String config, final int defaultValue)
    {
        return cmd != null ? Integer.parseInt(cmd.getOptionValue(config, String.valueOf(defaultValue))) : defaultValue;
    }

    public static void addCmdLineArgs(final Options options)
    {
        options.addOption(
                CFG_MIN_ALIGN_BASES, true, "Junction fragment min aligned bases, default: " + MIN_ALIGNMENT_BASES);

        options.addOption(CFG_MIN_JUNCTION_FRAGS, true, "Required fragments to keep a junction");

        options.addOption(CFG_MAX_JUNCTION_SUPPORTING_READS, true, "Limit to supporting reads added to a junction");
    }

    public ReadFilterConfig(
            final int minAlignmentBases, final int minMapQuality, final int minInsertAlignmentOverlap, final int minSoftClipLength,
            final int minSoftClipHighQual, final double minSoftClipHighQualPerc, final int minSupportingReadDistance,
            final int minIndelLength, final int minJunctionSupport, final int maxJunctionSupportingReads,
            final int maxDiscordantFragmentDistance)
    {
        MinAlignmentBases = minAlignmentBases;
        MinMapQuality = minMapQuality;
        MinInsertAlignmentOverlap = minInsertAlignmentOverlap;
        MinSoftClipLength = minSoftClipLength;
        MinIndelLength = minIndelLength;
        MinSoftClipHighQual = minSoftClipHighQual;
        MinSoftClipHighQualPerc = minSoftClipHighQualPerc;
        MinSupportingReadDistance = minSupportingReadDistance;
        MaxJunctionSupportingReads = maxJunctionSupportingReads;

        MinInsertLengthSupport = MIN_INSERT_LENGTH_SUPPORT;
        MinJunctionSupport = minJunctionSupport;
        MaxDiscordantFragmentDistance = maxDiscordantFragmentDistance;

        mFragmentLengthMin = DEFAULT_READ_LENGTH;
        mFragmentLengthMax = MAX_FRAGMENT_LENGTH;
    }

    public void setFragmentLengthMin(int minLength, int maxLength)
    {
        mFragmentLengthMin = minLength;
        mFragmentLengthMax = maxLength;
    }

    public int fragmentLengthMax() { return mFragmentLengthMax; }
}
