package com.hartwig.hmftools.purple.region;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.SegmentSupport;

import org.jetbrains.annotations.NotNull;

public final class ObservedRegionFile
{
    static final DecimalFormat FORMAT = new DecimalFormat("0.0000", new DecimalFormatSymbols(Locale.ENGLISH));

    private static final String EXTENSION = ".purple.observed.regions.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static void write(final String filename, final List<ObservedRegion> regions) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(regions));
    }

    @VisibleForTesting
    static List<String> toLines(final List<ObservedRegion> regions)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        regions.stream().map(ObservedRegionFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(TSV_DELIM)
                .add("chromosome")
                .add("start")
                .add("end")
                .add("ratioSupport")
                .add("support")
                .add("bafCount")
                .add("observedBAF")
                .add("depthWindowCount")
                .add("observedTumorRatio")
                .add("observedNormalRatio")
                .add("unnormalisedObservedNormalRatio")
                .add("germlineStatus")
                .add("svCluster")
                .add("gcContent")
                .add("minStart")
                .add("maxStart")
                .add("minorAlleleCopyNumberDeviation")
                .add("majorAlleleCopyNumberDeviation")
                .add("deviationPenalty")
                .add("eventPenalty")
                .add("refNormalisedCopyNumber")
                .add("tumorCopyNumber")
                .add("tumorBAF")
                .add("fittedTumorCopyNumber")
                .add("fittedBAF")
                .toString();
    }

    private static String toString(final ObservedRegion region)
    {
        return new StringJoiner(TSV_DELIM)
                .add(region.chromosome())
                .add(String.valueOf(region.start()))
                .add(String.valueOf(region.end()))
                .add(String.valueOf(region.ratioSupport()))
                .add(String.valueOf(region.support()))
                .add(String.valueOf(region.bafCount()))
                .add(FORMAT.format(region.observedBAF()))
                .add(String.valueOf(region.depthWindowCount()))
                .add(FORMAT.format(region.observedTumorRatio()))
                .add(FORMAT.format(region.observedNormalRatio()))
                .add(FORMAT.format(region.unnormalisedObservedNormalRatio()))
                .add(String.valueOf(region.germlineStatus()))
                .add(String.valueOf(region.svCluster()))
                .add(FORMAT.format(region.gcContent()))
                .add(String.valueOf(region.minStart()))
                .add(String.valueOf(region.maxStart()))
                .add(FORMAT.format(region.minorAlleleCopyNumberDeviation()))
                .add(FORMAT.format(region.majorAlleleCopyNumberDeviation()))
                .add(FORMAT.format(region.deviationPenalty()))
                .add(FORMAT.format(region.eventPenalty()))
                .add(FORMAT.format(region.refNormalisedCopyNumber()))
                .add(FORMAT.format(region.tumorCopyNumber()))
                .add(FORMAT.format(region.tumorBAF()))
                .add(FORMAT.format(region.fittedTumorCopyNumber()))
                .add(FORMAT.format(region.fittedBAF()))
                .toString();
    }

    public static List<ObservedRegion> read(final String fileName) throws IOException
    {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    @VisibleForTesting
    static List<ObservedRegion> fromLines(final List<String> lines)
    {
        List<ObservedRegion> regions = Lists.newArrayList();

        String header = lines.get(0);
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        lines.remove(0);

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            String chromosome = values[fieldsIndexMap.get("chromosome")];
            int start = Integer.parseInt(values[fieldsIndexMap.get("start")]);
            int end = Integer.parseInt(values[fieldsIndexMap.get("end")]);
            boolean ratioSupport = Boolean.parseBoolean(values[fieldsIndexMap.get("ratioSupport")]);
            SegmentSupport support = SegmentSupport.valueOf(values[fieldsIndexMap.get("support")]);
            int bafCount = Integer.parseInt(values[fieldsIndexMap.get("bafCount")]);
            double observedBAF = Double.parseDouble(values[fieldsIndexMap.get("observedBAF")]);
            int depthWindowCount = Integer.parseInt(values[fieldsIndexMap.get("depthWindowCount")]);
            double observedTumorRatio = Double.parseDouble(values[fieldsIndexMap.get("observedTumorRatio")]);
            double observedNormalRatio = Double.parseDouble(values[fieldsIndexMap.get("observedNormalRatio")]);
            double unnormalisedObservedNormalRatio = Double.parseDouble(values[fieldsIndexMap.get("unnormalisedObservedNormalRatio")]);
            GermlineStatus germlineStatus = GermlineStatus.valueOf(values[fieldsIndexMap.get("germlineStatus")]);
            boolean svCluster = Boolean.parseBoolean(values[fieldsIndexMap.get("svCluster")]);
            double gcContent = Double.parseDouble(values[fieldsIndexMap.get("gcContent")]);
            int minStart = Integer.parseInt(values[fieldsIndexMap.get("minStart")]);
            int maxStart = Integer.parseInt(values[fieldsIndexMap.get("maxStart")]);

            double minorAlleleCopyNumberDeviation = Double.parseDouble(values[fieldsIndexMap.get("minorAlleleCopyNumberDeviation")]);
            double majorAlleleCopyNumberDeviation = Double.parseDouble(values[fieldsIndexMap.get("majorAlleleCopyNumberDeviation")]);
            double deviationPenalty = Double.parseDouble(values[fieldsIndexMap.get("deviationPenalty")]);
            double eventPenalty = Double.parseDouble(values[fieldsIndexMap.get("eventPenalty")]);
            double refNormalisedCopyNumber = Double.parseDouble(values[fieldsIndexMap.get("refNormalisedCopyNumber")]);
            double tumorCopyNumber = Double.parseDouble(values[fieldsIndexMap.get("tumorCopyNumber")]);
            double tumorBAF = Double.parseDouble(values[fieldsIndexMap.get("tumorBAF")]);
            double fittedTumorCopyNumber = Double.parseDouble(values[fieldsIndexMap.get("fittedTumorCopyNumber")]);
            double fittedBAF = Double.parseDouble(values[fieldsIndexMap.get("fittedBAF")]);

            regions.add(new ObservedRegion(
                    chromosome, start, end, ratioSupport, support, bafCount, observedBAF, depthWindowCount,
                    observedTumorRatio, observedNormalRatio, unnormalisedObservedNormalRatio, germlineStatus, svCluster,
                    gcContent, minStart, maxStart, minorAlleleCopyNumberDeviation, majorAlleleCopyNumberDeviation,
                    deviationPenalty, eventPenalty, refNormalisedCopyNumber, tumorCopyNumber, tumorBAF,
                    fittedTumorCopyNumber, fittedBAF));
        }

        return regions;
    }
}
