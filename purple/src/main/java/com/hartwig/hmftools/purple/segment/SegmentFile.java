package com.hartwig.hmftools.purple.segment;

import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import org.jetbrains.annotations.NotNull;

public final class SegmentFile
{
    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");
    private static final String EXTENSION = ".purple.segment.tsv";
    private static final String DELIMITER = "\t";

    public static String generateFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static void write(final String filePath, Collection<ObservedRegion> fittedRegions) throws IOException
    {
        final Collection<String> lines = Lists.newArrayList();
        lines.add(header());
        fittedRegions.stream().map(SegmentFile::toString).forEach(lines::add);
        Files.write(new File(filePath).toPath(), lines);
    }

    public static List<ObservedRegion> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    private static List<ObservedRegion> fromLines(final List<String> lines)
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

        List<ObservedRegion> regions = Lists.newArrayList();

        for(final String line : lines)
        {
            String[] values = line.split(DELIMITER, -1);

            final ObservedRegion region = new ObservedRegion(
                    values[chrIndex],
                    Integer.parseInt(values[startIndex]),
                    Integer.parseInt(values[endIndex]),
                    Boolean.parseBoolean(values[ratioSupportIndex]),
                    SegmentSupport.valueOf(values[supportIndex]),
                    Integer.parseInt(values[bafCountIndex]),
                    Double.parseDouble(values[observedBAFIndex]),
                    Integer.parseInt(values[depthWindowCountIndex]),
                    Double.parseDouble(values[observedTumorRatioIndex]),
                    Double.parseDouble(values[observedNormalRatioIndex]),
                    Double.parseDouble(values[unnormalisedObservedNormalRatioIndex]),
                    GermlineStatus.valueOf(values[germlineStatusIndex]),
                    false,
                    Double.parseDouble(values[gcContentIndex]),
                    Integer.parseInt(values[minStartIndex]),
                    Integer.parseInt(values[maxStartIndex]),
                    Double.parseDouble(values[minorAlleleCopyNumberDeviationIndex]),
                    Double.parseDouble(values[majorAlleleCopyNumberDeviationIndex]),
                    Double.parseDouble(values[deviationPenaltyIndex]),
                    Double.parseDouble(values[eventPenaltyIndex]),
                    Double.parseDouble(values[refNormalisedCopyNumberIndex]),
                    Double.parseDouble(values[tumorCopyNumberIndex]),
                    Double.parseDouble(values[tumorBAFIndex]),
                    Double.parseDouble(values[fittedTumorCopyNumberIndex]),
                    Double.parseDouble(values[fittedBAFIndex]));

            regions.add(region);
        }

        return regions;
    }

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

    private static String toString(final ObservedRegion region)
    {
        return new StringJoiner(DELIMITER)
                .add(region.chromosome())
                .add(String.valueOf(region.start()))
                .add(String.valueOf(region.end()))
                .add(String.valueOf(region.germlineStatus()))
                .add(String.valueOf(region.bafCount()))
                .add(FORMAT.format(region.observedBAF()))
                .add(FORMAT.format(region.minorAlleleCopyNumber()))
                .add(FORMAT.format(region.minorAlleleCopyNumberDeviation()))
                .add(FORMAT.format(region.observedTumorRatio()))
                .add(FORMAT.format(region.observedNormalRatio()))
                .add(FORMAT.format(region.unnormalisedObservedNormalRatio()))
                .add(FORMAT.format(region.majorAlleleCopyNumber()))
                .add(FORMAT.format(region.majorAlleleCopyNumberDeviation()))
                .add(FORMAT.format(region.deviationPenalty()))
                .add(FORMAT.format(region.tumorCopyNumber()))
                .add(FORMAT.format(region.fittedTumorCopyNumber()))
                .add(FORMAT.format(region.fittedBAF()))
                .add(FORMAT.format(region.refNormalisedCopyNumber()))
                .add(String.valueOf(region.ratioSupport()))
                .add(String.valueOf(region.support()))
                .add(String.valueOf(region.depthWindowCount()))
                .add(FORMAT.format(region.tumorBAF()))
                .add(FORMAT.format(region.gcContent()))
                .add(FORMAT.format(region.eventPenalty()))
                .add(String.valueOf(region.minStart()))
                .add(String.valueOf(region.maxStart()))
                .toString();
    }
}
