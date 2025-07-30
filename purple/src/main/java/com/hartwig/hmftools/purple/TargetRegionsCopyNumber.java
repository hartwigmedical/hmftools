package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.TaggedRegion;

public record TargetRegionsCopyNumber(CobaltRatio mCobaltRatio, List<TaggedRegion> mOverlappingRegions)
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
                .toString();
    }

    public String toTSV()
    {
        BaseRegion cobaltRegion = mCobaltRatio.baseRegion();
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
                .toString();
    }
}
