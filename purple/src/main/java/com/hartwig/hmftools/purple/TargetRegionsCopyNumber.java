package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.TaggedRegion;

public record TargetRegionsCopyNumber(CobaltRatio mCobaltRatio,
                                      List<TaggedRegion> mOverlappingRegions,
                                      PurpleCopyNumber mPurpleCopyNumber)
{
    public static String tsvFileHeader()
    {
        return new StringJoiner("\t")
                .add("chromosome")
                .add("windowStart")
                .add("windowEnd")
                .add("panelRegions")
                .add("masked")
                .add("averageDepth")
                .add("windowGCContent")
                .add("windowTumorRatio")
                .add("regionStart")
                .add("regionEnd")
                .add("copyNumber")
                .add("minorAlleleCopNumber")
                .add("depthWindowCount")
                .add("bafCount")
                .add("GCContent")
                .toString();
    }

    public String toTSV()
    {
        BaseRegion cobaltRegion = mCobaltRatio.window();
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
                .add(String.valueOf(mCobaltRatio.tumorReadDepth()))
                .add(String.valueOf(mCobaltRatio.tumorGcContent()))
                .add(String.valueOf(mCobaltRatio.tumorGCRatio()))
                .add(String.valueOf(mPurpleCopyNumber.start()))
                .add(String.valueOf(mPurpleCopyNumber.end()))
                .add(String.valueOf(mPurpleCopyNumber.averageTumorCopyNumber()))
                .add(String.valueOf(mPurpleCopyNumber.minorAlleleCopyNumber()))
                .add(String.valueOf(mPurpleCopyNumber.depthWindowCount()))
                .add(String.valueOf(mPurpleCopyNumber.bafCount()))
                .add(String.valueOf(mPurpleCopyNumber.gcContent()))
                .toString();
    }
}
