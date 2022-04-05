package com.hartwig.hmftools.purple.segment;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.purple.region.ImmutableFittedRegion;

import org.jetbrains.annotations.NotNull;

public final class SegmentFile
{
    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");
    private static final String EXTENSION = ".purple.segment.tsv";
    private static final String DELIMITER = "\t";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static void write(@NotNull final String filePath, @NotNull Collection<FittedRegion> fittedRegions) throws IOException
    {
        final Collection<String> lines = Lists.newArrayList();
        lines.add(header());
        fittedRegions.stream().map(SegmentFile::toString).forEach(lines::add);
        Files.write(new File(filePath).toPath(), lines);
    }

    @NotNull
    public static List<FittedRegion> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    private static List<FittedRegion> fromLines(final List<String> lines)
    {
        final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
        lines.remove(0);

        int chrIndex = fieldsIndexMap.get("chromosome");
        int startIndex = fieldsIndexMap.get("start");
        int endIndex = fieldsIndexMap.get("end");
        int germlineStatusIndex = fieldsIndexMap.get("germlineStatus");
        int bafCountIndex = fieldsIndexMap.get("bafCount");
        int observedBAFIndex = fieldsIndexMap.get("observedBAF");
        int minorAlleleCopyNumberDeviationIndex = fieldsIndexMap.get("minorAlleleCopyNumberDeviation");
        int observedTumorRatioIndex = fieldsIndexMap.get("observedTumorRatio");
        int observedNormalRatioIndex = fieldsIndexMap.get("observedNormalRatio");
        int unnormalisedObservedNormalRatioIndex = fieldsIndexMap.get("unnormalisedObservedNormalRatio");
        int majorAlleleCopyNumberDeviationIndex = fieldsIndexMap.get("majorAlleleCopyNumberDeviation");
        int deviationPenaltyIndex = fieldsIndexMap.get("deviationPenalty");
        int tumorCopyNumberIndex = fieldsIndexMap.get("tumorCopyNumber");
        int fittedTumorCopyNumberIndex = fieldsIndexMap.get("fittedTumorCopyNumber");
        int fittedBAFIndex = fieldsIndexMap.get("fittedBAF");
        int refNormalisedCopyNumberIndex = fieldsIndexMap.get("refNormalisedCopyNumber");
        int ratioSupportIndex = fieldsIndexMap.get("ratioSupport");
        int supportIndex = fieldsIndexMap.get("support");
        int depthWindowCountIndex = fieldsIndexMap.get("depthWindowCount");
        int tumorBAFIndex = fieldsIndexMap.get("tumorBAF");
        int gcContentIndex = fieldsIndexMap.get("gcContent");
        int eventPenaltyIndex = fieldsIndexMap.get("eventPenalty");
        int minStartIndex = fieldsIndexMap.get("minStart");
        int maxStartIndex = fieldsIndexMap.get("maxStart");

        List<FittedRegion> regions = Lists.newArrayList();

        for(final String line : lines)
        {
            String[] values = line.split(DELIMITER, -1);

            final ImmutableFittedRegion.Builder builder = ImmutableFittedRegion.builder()
                    .chromosome(values[chrIndex])
                    .start(Integer.parseInt(values[startIndex]))
                    .end(Integer.parseInt(values[endIndex]))
                    .germlineStatus(GermlineStatus.valueOf(values[germlineStatusIndex]))
                    .bafCount(Integer.parseInt(values[bafCountIndex]))
                    .observedBAF(Double.parseDouble(values[observedBAFIndex]))
                    .minorAlleleCopyNumberDeviation(Double.parseDouble(values[minorAlleleCopyNumberDeviationIndex]))
                    .observedTumorRatio(Double.parseDouble(values[observedTumorRatioIndex]))
                    .observedNormalRatio(Double.parseDouble(values[observedNormalRatioIndex]))
                    .unnormalisedObservedNormalRatio(Double.parseDouble(values[unnormalisedObservedNormalRatioIndex]))
                    .majorAlleleCopyNumberDeviation(Double.parseDouble(values[majorAlleleCopyNumberDeviationIndex]))
                    .deviationPenalty(Double.parseDouble(values[deviationPenaltyIndex]))
                    .tumorCopyNumber(Double.parseDouble(values[tumorCopyNumberIndex]))
                    .fittedTumorCopyNumber(Double.parseDouble(values[fittedTumorCopyNumberIndex]))
                    .fittedBAF(Double.parseDouble(values[fittedBAFIndex]))
                    .refNormalisedCopyNumber(Double.parseDouble(values[refNormalisedCopyNumberIndex]))
                    .ratioSupport(Boolean.parseBoolean(values[ratioSupportIndex]))
                    .support(SegmentSupport.valueOf(values[supportIndex]))
                    .depthWindowCount(Integer.parseInt(values[depthWindowCountIndex]))
                    .tumorBAF(Double.parseDouble(values[tumorBAFIndex]))
                    .gcContent(Double.parseDouble(values[gcContentIndex]))
                    .eventPenalty(Double.parseDouble(values[eventPenaltyIndex]))
                    .minStart(Integer.parseInt(values[minStartIndex]))
                    .maxStart(Integer.parseInt(values[maxStartIndex]))
                    .svCluster(false);

            regions.add(builder.build());
        }

        return regions;
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER, "", "")
                .add("chromosome")
                .add("start")
                .add("end")
                .add("germlineStatus")
                .add("bafCount")
                .add("observedBAF")
                .add("minorAlleleCopyNumber")
                .add("minorAlleleCopyNumberDeviation")
                .add("observedTumorRatio")
                .add("observedNormalRatio")
                .add("unnormalisedObservedNormalRatio")
                .add("majorAlleleCopyNumber")
                .add("majorAlleleCopyNumberDeviation")
                .add("deviationPenalty")
                .add("tumorCopyNumber")
                .add("fittedTumorCopyNumber")
                .add("fittedBAF")
                .add("refNormalisedCopyNumber")
                .add("ratioSupport")
                .add("support")
                .add("depthWindowCount")
                .add("tumorBAF")
                .add("gcContent")
                .add("eventPenalty")
                .add("minStart")
                .add("maxStart")
                .toString();
    }

    @NotNull
    static String toString(@NotNull final FittedRegion copyNumber)
    {
        return new StringJoiner(DELIMITER)
                .add(copyNumber.chromosome())
                .add(String.valueOf(copyNumber.start()))
                .add(String.valueOf(copyNumber.end()))
                .add(String.valueOf(copyNumber.germlineStatus()))
                .add(String.valueOf(copyNumber.bafCount()))
                .add(FORMAT.format(copyNumber.observedBAF()))
                .add(FORMAT.format(copyNumber.minorAlleleCopyNumber()))
                .add(FORMAT.format(copyNumber.minorAlleleCopyNumberDeviation()))
                .add(FORMAT.format(copyNumber.observedTumorRatio()))
                .add(FORMAT.format(copyNumber.observedNormalRatio()))
                .add(FORMAT.format(copyNumber.unnormalisedObservedNormalRatio()))
                .add(FORMAT.format(copyNumber.majorAlleleCopyNumber()))
                .add(FORMAT.format(copyNumber.majorAlleleCopyNumberDeviation()))
                .add(FORMAT.format(copyNumber.deviationPenalty()))
                .add(FORMAT.format(copyNumber.tumorCopyNumber()))
                .add(FORMAT.format(copyNumber.fittedTumorCopyNumber()))
                .add(FORMAT.format(copyNumber.fittedBAF()))
                .add(FORMAT.format(copyNumber.refNormalisedCopyNumber()))
                .add(String.valueOf(copyNumber.ratioSupport()))
                .add(String.valueOf(copyNumber.support()))
                .add(String.valueOf(copyNumber.depthWindowCount()))
                .add(FORMAT.format(copyNumber.tumorBAF()))
                .add(FORMAT.format(copyNumber.gcContent()))
                .add(FORMAT.format(copyNumber.eventPenalty()))
                .add(String.valueOf(copyNumber.minStart()))
                .add(String.valueOf(copyNumber.maxStart()))
                .toString();
    }
}
