package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

public final class PurpleSegment
{
    public final String Chromosome;
    public final int PosStart;
    public final int PosEnd;

    public final boolean RatioSupport;

    public final SegmentSupport Support;

    public final int BafCount;
    public final double ObservedBAF;
    public final int DepthWindowCount;
    public final double ObservedTumorRatio;
    public final double ObservedNormalRatio;
    public final double UnnormalisedObservedNormalRatio;

    public final GermlineStatus GermlineState;
    public final boolean SvCluster;
    public final double GcContent;

    public final int MinStart;
    public final int MaxStart;

    // fitted region fields
    public final double MinorAlleleCopyNumberDeviation;
    public final double MajorAlleleCopyNumberDeviation;
    public final double DeviationPenalty;
    public final double EventPenalty;
    public final double RefNormalisedCopyNumber;
    public final double TumorCopyNumber;
    public final double TumorBAF;
    public final double FittedTumorCopyNumber;
    public final double FittedBAF;

    public static final String EXTENSION = ".purple.segment.tsv";

    public PurpleSegment(
            final String chromosome, final int posStart, final int posEnd, final boolean ratioSupport,
            final SegmentSupport support, final int bafCount, final double observedBAF, final int depthWindowCount,
            final double observedTumorRatio, final double observedNormalRatio, final double unnormalisedObservedNormalRatio,
            final GermlineStatus germlineStatus, final boolean svCluster, final double gcContent, final int minStart, final int maxStart,
            double minorAlleleCopyNumberDeviation, double majorAlleleCopyNumberDeviation, double deviationPenalty, double eventPenalty,
            double refNormalisedCopyNumber, double tumorCopyNumber, double tumorBAF, double fittedTumorCopyNumber, double fittedBAF)
    {
        Chromosome = chromosome;
        PosStart = posStart;
        PosEnd = posEnd;

        RatioSupport = ratioSupport;
        Support = support;
        BafCount = bafCount;
        ObservedBAF = observedBAF;
        DepthWindowCount = depthWindowCount;
        ObservedTumorRatio = observedTumorRatio;
        ObservedNormalRatio = observedNormalRatio;
        UnnormalisedObservedNormalRatio = unnormalisedObservedNormalRatio;
        GermlineState = germlineStatus;
        SvCluster = svCluster;
        GcContent = gcContent;
        MinStart = minStart;
        MaxStart = maxStart;

        MinorAlleleCopyNumberDeviation = minorAlleleCopyNumberDeviation;
        MajorAlleleCopyNumberDeviation = majorAlleleCopyNumberDeviation;
        DeviationPenalty = deviationPenalty;
        EventPenalty = eventPenalty;
        RefNormalisedCopyNumber = refNormalisedCopyNumber;
        TumorCopyNumber = tumorCopyNumber;
        TumorBAF = tumorBAF;
        FittedTumorCopyNumber = fittedTumorCopyNumber;
        FittedBAF = fittedBAF;
    }

    public double minorAlleleCopyNumber() { return TumorCopyNumber - majorAlleleCopyNumber(); }
    public double majorAlleleCopyNumber() { return TumorBAF * TumorCopyNumber; }

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + EXTENSION;
    }

    public static List<PurpleSegment> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    private static List<PurpleSegment> fromLines(final List<String> lines)
    {
        final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);
        lines.remove(0);

        int chrIndex = fieldsIndexMap.get("chromosome");
        int startIndex = fieldsIndexMap.get("start");
        int endIndex = fieldsIndexMap.get("end");
        int germlineStatusIndex = fieldsIndexMap.get("germlineStatus");
        int svClusterIndex = fieldsIndexMap.get("svCluster");
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

        List<PurpleSegment> regions = Lists.newArrayList();

        for(final String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            PurpleSegment region = new PurpleSegment(
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
                    Boolean.parseBoolean(values[svClusterIndex]),
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

    public static void write(final String filename, final List<PurpleSegment> segments) throws IOException
    {
        final Collection<String> lines = Lists.newArrayList();
        lines.add(PurpleSegment.header());

        segments.stream().map(x -> x.toTsv()).forEach(lines::add);
        Files.write(new File(filename).toPath(), lines);
    }

    private static String header()
    {
        return new StringJoiner(TSV_DELIM, "", "")
                .add("chromosome")
                .add("start")
                .add("end")
                .add("germlineStatus")
                .add("svCluster")
                .add("bafCount")
                .add("observedBAF")
                .add("minorAlleleCopyNumberDeviation")
                .add("observedTumorRatio")
                .add("observedNormalRatio")
                .add("unnormalisedObservedNormalRatio")
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

    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");

    private String toTsv()
    {
        return new StringJoiner(TSV_DELIM)
                .add(Chromosome)
                .add(String.valueOf(PosStart))
                .add(String.valueOf(PosEnd))
                .add(String.valueOf(GermlineState))
                .add(String.valueOf(SvCluster))
                .add(String.valueOf(BafCount))
                .add(FORMAT.format(ObservedBAF))
                .add(FORMAT.format(MinorAlleleCopyNumberDeviation))
                .add(FORMAT.format(ObservedTumorRatio))
                .add(FORMAT.format(ObservedNormalRatio))
                .add(FORMAT.format(UnnormalisedObservedNormalRatio))
                .add(FORMAT.format(MajorAlleleCopyNumberDeviation))
                .add(FORMAT.format(DeviationPenalty))
                .add(FORMAT.format(TumorCopyNumber))
                .add(FORMAT.format(FittedTumorCopyNumber))
                .add(FORMAT.format(FittedBAF))
                .add(FORMAT.format(RefNormalisedCopyNumber))
                .add(String.valueOf(RatioSupport))
                .add(String.valueOf(Support))
                .add(String.valueOf(DepthWindowCount))
                .add(FORMAT.format(TumorBAF))
                .add(FORMAT.format(GcContent))
                .add(FORMAT.format(EventPenalty))
                .add(String.valueOf(MinStart))
                .add(String.valueOf(MaxStart))
                .toString();
    }
}
