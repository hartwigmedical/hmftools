package com.hartwig.hmftools.purple.targeted;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Collection;
import java.util.List;
import java.util.Locale;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.TaggedRegion;

public class TargetRegionsCopyNumberFile
{
    public static final String EXTENSION = ".purple.target_region_cn.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + EXTENSION;
    }

    public static void write(final String filename, final List<TargetRegionsCopyNumber> segments) throws IOException
    {
        final Collection<String> lines = Lists.newArrayList();
        lines.add(tsvFileHeader());
        segments.forEach(x -> lines.add(toTSV(x)));
        Files.write(new File(filename).toPath(), lines);
    }

    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000", new DecimalFormatSymbols(Locale.ENGLISH));

    private static String tsvFileHeader()
    {
        return new StringJoiner(TSV_DELIM)
                .add("chromosome").add("windowStart").add("windowEnd").add("bedRegions")
                .add("masked").add("averageDepth").add("windowGCContent").add("windowTumorRatio")
                .add("regionStart").add("regionEnd").add("copyNumber").add("minorAlleleCopyNumber").add("depthWindowCount")
                .add("bafCount").add("germlineStatus").add("copyNumberMethod").toString();
    }

    private static String toTSV(final TargetRegionsCopyNumber data)
    {
        CobaltRatio cobaltRatio = data.cobaltRatio();
        ChrBaseRegion cobaltRegion = cobaltRatio.window();
        StringJoiner panelRegionsStringJoiner = new StringJoiner(ITEM_DELIM);

        for(TaggedRegion region : data.overlappingRegions())
        {
            panelRegionsStringJoiner.add(region.formatted());
        }

        boolean masked = cobaltRatio.tumorGCRatio() < 0;

        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(cobaltRatio.chromosome()).add(String.valueOf(cobaltRegion.start())).add(String.valueOf(cobaltRegion.end()));
        sj.add(panelRegionsStringJoiner.toString()).add(String.valueOf(masked));

        sj.add(FORMAT.format(cobaltRatio.tumorReadDepth()));
        sj.add(FORMAT.format(cobaltRatio.tumorGcContent()));
        sj.add(FORMAT.format(cobaltRatio.tumorGCRatio()));

        PurpleCopyNumber copyNumber = data.purpleCopyNumber();
        sj.add(String.valueOf(copyNumber.start()));
        sj.add(String.valueOf(copyNumber.end()));
        sj.add(FORMAT.format(copyNumber.averageTumorCopyNumber()));
        sj.add(FORMAT.format(copyNumber.minorAlleleCopyNumber()));
        sj.add(String.valueOf(copyNumber.depthWindowCount()));
        sj.add(String.valueOf(copyNumber.bafCount()));
        sj.add(data.germlineStatus().name());
        sj.add(copyNumber.method().name());

        return sj.toString();
    }
}
