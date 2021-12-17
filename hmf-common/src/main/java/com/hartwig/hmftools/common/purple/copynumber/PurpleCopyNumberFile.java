package com.hartwig.hmftools.common.purple.copynumber;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;

public final class PurpleCopyNumberFile
{
    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");
    private static final String DELIMITER = "\t";

    private static final String SOMATIC_EXTENSION = ".purple.cnv.somatic.tsv";
    private static final String SOMATIC_EXTENSION_OLD = ".purple.cnv";

    private PurpleCopyNumberFile()
    {
    }

    @NotNull
    public static String generateFilenameForWriting(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + SOMATIC_EXTENSION;
    }

    @NotNull
    public static String generateFilenameForReading(@NotNull final String basePath, @NotNull final String sample)
    {
        String filename = basePath + File.separator + sample + SOMATIC_EXTENSION;
        return (new File(filename).exists()) ? filename : basePath + File.separator + sample + SOMATIC_EXTENSION_OLD;
    }

    @NotNull
    public static List<PurpleCopyNumber> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<PurpleCopyNumber> copyNumbers) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(copyNumbers));
    }

    @VisibleForTesting
    @NotNull
    public static List<String> toLines(@NotNull final List<PurpleCopyNumber> copyNumbers)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        copyNumbers.stream().map(PurpleCopyNumberFile::toString).forEach(lines::add);
        return lines;
    }

    @VisibleForTesting
    @NotNull
    public static List<PurpleCopyNumber> fromLines(final List<String> lines)
    {
        final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
        lines.remove(0);

        int chrIndex = fieldsIndexMap.get("chromosome");
        int startIndex = fieldsIndexMap.get("start");
        int endIndex = fieldsIndexMap.get("end");
        int cnIndex = fieldsIndexMap.get("copyNumber");
        int bafCountIndex = fieldsIndexMap.get("bafCount");
        int observedBAFIndex = fieldsIndexMap.get("observedBAF");
        int bafIndex = fieldsIndexMap.get("baf");
        int segmentStartSupportIndex = fieldsIndexMap.get("segmentStartSupport");
        int segmentEndSupportIndex = fieldsIndexMap.get("segmentEndSupport");
        int methodIndex = fieldsIndexMap.get("method");
        int depthWindowCountIndex = fieldsIndexMap.get("depthWindowCount");
        int gcContentIndex = fieldsIndexMap.get("gcContent");
        int minStartIndex = fieldsIndexMap.get("minStart");
        int maxStartIndex = fieldsIndexMap.get("maxStart");

        List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();

        for(final String line : lines)
        {
            String[] values = line.split(DELIMITER, -1);

            final ImmutablePurpleCopyNumber.Builder builder = ImmutablePurpleCopyNumber.builder()
                    .chromosome(values[chrIndex])
                    .start(Integer.parseInt(values[startIndex]))
                    .end(Integer.parseInt(values[endIndex]))
                    .bafCount(Integer.parseInt(values[bafCountIndex]))
                    .averageTumorCopyNumber(Double.parseDouble(values[cnIndex]))
                    .averageObservedBAF(Double.parseDouble(values[observedBAFIndex]))
                    .averageActualBAF(Double.parseDouble(values[bafIndex]))
                    .segmentStartSupport(SegmentSupport.valueOf(values[segmentStartSupportIndex]))
                    .segmentEndSupport(SegmentSupport.valueOf(values[segmentEndSupportIndex]))
                    .method(CopyNumberMethod.valueOf(values[methodIndex]))
                    .gcContent(Double.parseDouble(values[gcContentIndex]))
                    .depthWindowCount(Integer.parseInt(values[depthWindowCountIndex]))
                    .minStart(Long.parseLong(values[minStartIndex]))
                    .maxStart(Long.parseLong(values[maxStartIndex]));

            copyNumbers.add(builder.build());
        }

        return copyNumbers;
    }

    private static String header()
    {
        return new StringJoiner(DELIMITER, "", "")
                .add("chromosome")
                .add("start")
                .add("end")
                .add("copyNumber")
                .add("bafCount")
                .add("observedBAF")
                .add("baf")
                .add("segmentStartSupport")
                .add("segmentEndSupport")
                .add("method")
                .add("depthWindowCount")
                .add("gcContent")
                .add("minStart")
                .add("maxStart")
                .add("minorAlleleCopyNumber")
                .add("majorAlleleCopyNumber")
                .toString();
    }

    private static String toString(@NotNull final PurpleCopyNumber copyNumber)
    {
        return new StringJoiner(DELIMITER).add(copyNumber.chromosome())
                .add(String.valueOf(copyNumber.start()))
                .add(String.valueOf(copyNumber.end()))
                .add(FORMAT.format(copyNumber.averageTumorCopyNumber()))
                .add(String.valueOf(copyNumber.bafCount()))
                .add(FORMAT.format(copyNumber.averageObservedBAF()))
                .add(FORMAT.format(copyNumber.averageActualBAF()))
                .add(String.valueOf(copyNumber.segmentStartSupport()))
                .add(String.valueOf(copyNumber.segmentEndSupport()))
                .add(String.valueOf(copyNumber.method()))
                .add(String.valueOf(copyNumber.depthWindowCount()))
                .add(FORMAT.format(copyNumber.gcContent()))
                .add(String.valueOf(copyNumber.minStart()))
                .add(String.valueOf(copyNumber.maxStart()))
                .add(FORMAT.format(copyNumber.minorAlleleCopyNumber()))
                .add(FORMAT.format(copyNumber.majorAlleleCopyNumber()))
                .toString();
    }

    private static PurpleCopyNumber fromString(@NotNull final String copyNumber)
    {
        String[] values = copyNumber.split(DELIMITER);
        final ImmutablePurpleCopyNumber.Builder builder = ImmutablePurpleCopyNumber.builder()
                .chromosome(values[0])
                .start(Integer.parseInt(values[1]))
                .end(Integer.parseInt(values[2]))
                .averageTumorCopyNumber(Double.parseDouble(values[3]))
                .bafCount(Integer.parseInt(values[4]))
                .averageObservedBAF(Double.parseDouble(values[5]))
                .averageActualBAF(Double.parseDouble(values[6]))
                .segmentStartSupport(SegmentSupport.valueOf(values[7]))
                .segmentEndSupport(SegmentSupport.valueOf(values[8]))
                .method(CopyNumberMethod.valueOf(values[9]))
                .depthWindowCount(Integer.parseInt(values[10]))
                .gcContent(Double.parseDouble(values[11]))
                .minStart(Long.parseLong(values[12]))
                .maxStart(Long.parseLong(values[13]));

        return builder.build();
    }
}
