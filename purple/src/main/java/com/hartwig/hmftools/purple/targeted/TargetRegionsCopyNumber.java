package com.hartwig.hmftools.purple.targeted;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.List;
import java.util.Locale;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.TaggedRegion;

public record TargetRegionsCopyNumber(
        CobaltRatio cobaltRatio, List<TaggedRegion> overlappingRegions, PurpleCopyNumber purpleCopyNumber, GermlineStatus germlineStatus)
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
                .add("germlineStatus")
                .add("copyNumberMethod")
                .toString();
    }

    public String toTSV()
    {
        ChrBaseRegion cobaltRegion = cobaltRatio.window();
        StringJoiner panelRegionsStringJoiner = new StringJoiner(ITEM_DELIM);

        for(TaggedRegion region : overlappingRegions)
        {
            panelRegionsStringJoiner.add(region.formatted());
        }

        boolean masked = cobaltRatio.tumorGCRatio() < 0;
        return new StringJoiner(TSV_DELIM).add(cobaltRatio.chromosome())
                .add(String.valueOf(cobaltRegion.start()))
                .add(String.valueOf(cobaltRegion.end()))
                .add(panelRegionsStringJoiner.toString())
                .add(String.valueOf(masked))
                .add(FORMAT.format(cobaltRatio.tumorReadDepth()))
                .add(FORMAT.format(cobaltRatio.tumorGcContent()))
                .add(FORMAT.format(cobaltRatio.tumorGCRatio()))
                .add(String.valueOf(purpleCopyNumber.start()))
                .add(String.valueOf(purpleCopyNumber.end()))
                .add(FORMAT.format(purpleCopyNumber.averageTumorCopyNumber()))
                .add(FORMAT.format(purpleCopyNumber.minorAlleleCopyNumber()))
                .add(String.valueOf(purpleCopyNumber.depthWindowCount()))
                .add(String.valueOf(purpleCopyNumber.bafCount()))
                .add(germlineStatus.name())
                .add(purpleCopyNumber.method().name())
                .toString();
    }
}
