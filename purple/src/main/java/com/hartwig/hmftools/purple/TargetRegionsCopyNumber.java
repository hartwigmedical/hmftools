package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.List;
import java.util.Locale;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.TaggedRegion;

public record TargetRegionsCopyNumber(CobaltRatio mCobaltRatio,
                                      List<TaggedRegion> mOverlappingRegions,
                                      PurpleCopyNumber mPurpleCopyNumber)
{
    static final DecimalFormat FORMAT = new DecimalFormat("0.0000", new DecimalFormatSymbols(Locale.ENGLISH));

    public static String tsvFileHeader()
    {
        return new StringJoiner("\t")
                .add("chromosome")
                .add("windowStart")
                .add("windowEnd")
                .add("bedRegions")
                .add("masked")
                .add("averageDepth")
                .add("windowGCContent")
                .add("windowTumorRatio")
                .add("regionStart")
                .add("regionEnd")
                .add("copyNumber")
                .add("minorAlleleCopyNumber")
                .add("depthWindowCount")
                .add("bafCount")
                .add("GCContent")
                .add("CNMethod")
                .toString();
    }

    public String toTSV()
    {
        ChrBaseRegion cobaltRegion = mCobaltRatio.window();
        StringJoiner panelRegionsStringJoiner = new StringJoiner(":");
        for(TaggedRegion region : mOverlappingRegions)
        {
            panelRegionsStringJoiner.add(region.formatted());
        }
        boolean masked = mCobaltRatio.tumorGCRatio() < 0;
        return new StringJoiner(TSV_DELIM).add(mCobaltRatio.chromosome())
                .add(String.valueOf(cobaltRegion.start()))
                .add(String.valueOf(cobaltRegion.end()))
                .add(panelRegionsStringJoiner.toString())
                .add(String.valueOf(masked))
                .add(FORMAT.format(mCobaltRatio.tumorReadDepth()))
                .add(FORMAT.format(mCobaltRatio.tumorGcContent()))
                .add(FORMAT.format(mCobaltRatio.tumorGCRatio()))
                .add(String.valueOf(mPurpleCopyNumber.start()))
                .add(String.valueOf(mPurpleCopyNumber.end()))
                .add(FORMAT.format(mPurpleCopyNumber.averageTumorCopyNumber()))
                .add(FORMAT.format(mPurpleCopyNumber.minorAlleleCopyNumber()))
                .add(String.valueOf(mPurpleCopyNumber.depthWindowCount()))
                .add(String.valueOf(mPurpleCopyNumber.bafCount()))
                .add(FORMAT.format(mPurpleCopyNumber.gcContent()))
                .add(mPurpleCopyNumber.method().name())
                .toString();
    }
}
