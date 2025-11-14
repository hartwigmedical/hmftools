package com.hartwig.hmftools.common.linx;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.TranscriptCodingType;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;
import com.hartwig.hmftools.common.purple.ReportedStatus;

import org.immutables.value.Value;

@Value.Immutable
public abstract class LinxBreakend
{
    public abstract int id();
    public abstract int svId();
    public abstract String vcfId();
    public abstract String coords();
    public abstract boolean isStart();
    public abstract String gene();
    public abstract String transcriptId();
    public abstract boolean canonical();
    public abstract String geneOrientation();
    public abstract boolean disruptive();
    public abstract ReportedStatus reportedStatus();
    public abstract double undisruptedCopyNumber();
    public abstract TranscriptRegionType regionType();
    public abstract TranscriptCodingType codingType();
    public abstract String biotype();
    public abstract int exonUp();
    public abstract int exonDown();
    public abstract int exonicBasePhase();
    public abstract int nextSpliceExonRank();
    public abstract int nextSpliceExonPhase();
    public abstract int nextSpliceDistance();
    public abstract int totalExonCount();

    private static final String FILE_EXTENSION = ".linx.breakend.tsv";
    private static final String GERMLINE_FILE_EXTENSION = ".linx.germline.breakend.tsv";

    public static final String BREAKEND_ORIENTATION_UPSTREAM = "Upstream";
    public static final String BREAKEND_ORIENTATION_DOWNSTREAM = "Downstream";

    public static final String BREAKEND_COORD_DELIM = ":";

    public static String generateFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static String generateFilename(final String basePath, final String sample, boolean isGermline)
    {
        return basePath + File.separator + sample + (isGermline ? GERMLINE_FILE_EXTENSION : FILE_EXTENSION);
    }

    public static List<LinxBreakend> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(final String filename, List<LinxBreakend> breakends) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(breakends));
    }

    public static String coordsStr(final String chromosome, final int position, final byte orientation)
    {
        return new StringJoiner(BREAKEND_COORD_DELIM)
                .add(chromosome)
                .add(String.valueOf(position))
                .add(String.valueOf(orientation)).toString();
    }

    public static String chromosomeFromCoords(final String coordStr)
    {
        return coordStr.split(BREAKEND_COORD_DELIM, 3)[0];
    }

    public static Integer positionFromCoords(final String coordStr)
    {
        return Integer.parseInt(coordStr.split(BREAKEND_COORD_DELIM, 3)[1]);
    }

    public static byte orientationFromCoords(final String coordStr)
    {
        return Byte.parseByte(coordStr.split(BREAKEND_COORD_DELIM, 3)[2]);
    }

    private static List<String> toLines(final List<LinxBreakend> breakends)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        breakends.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    private static List<LinxBreakend> fromLines(final List<String> lines)
    {
        String header = lines.get(0);
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        lines.remove(0);

        List<LinxBreakend> breakends = Lists.newArrayList();

        Integer vcfIdIndex = fieldsIndexMap.get("vcfId");
        Integer coordsIndex = fieldsIndexMap.get("coords");
        Integer reportedIndex = fieldsIndexMap.get("reportedDisruption");
        Integer reportedStatusIndex = fieldsIndexMap.get("reportedStatus");

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM);

            ReportedStatus reportedStatus;

            if(reportedStatusIndex != null)
                reportedStatus = ReportedStatus.valueOf(values[reportedStatusIndex]);
            else
                reportedStatus = Boolean.parseBoolean(values[reportedIndex]) ? ReportedStatus.REPORTED : ReportedStatus.NONE;

            breakends.add(ImmutableLinxBreakend.builder()
                    .id(Integer.parseInt(values[fieldsIndexMap.get("id")]))
                    .svId(Integer.parseInt(values[fieldsIndexMap.get("svId")]))
                    .vcfId(vcfIdIndex != null ? values[vcfIdIndex] : "")
                    .coords(coordsIndex != null ? values[coordsIndex] : "")
                    .isStart(Boolean.parseBoolean(values[fieldsIndexMap.get("isStart")]))
                    .gene(values[fieldsIndexMap.get("gene")])
                    .transcriptId(values[fieldsIndexMap.get("transcriptId")])
                    .canonical(Boolean.parseBoolean(values[fieldsIndexMap.get("canonical")]))
                    .geneOrientation(values[fieldsIndexMap.get("geneOrientation")])
                    .disruptive(Boolean.parseBoolean(values[fieldsIndexMap.get("disruptive")]))
                    .reportedStatus(reportedStatus)
                    .undisruptedCopyNumber(Double.parseDouble(values[fieldsIndexMap.get("undisruptedCopyNumber")]))
                    .regionType(TranscriptRegionType.valueOf(values[fieldsIndexMap.get("regionType")]))
                    .codingType(TranscriptCodingType.valueOf(values[fieldsIndexMap.get("codingType")]))
                    .biotype(values[fieldsIndexMap.get("biotype")])
                    .exonicBasePhase(Integer.parseInt(values[fieldsIndexMap.get("exonicBasePhase")]))
                    .nextSpliceExonRank(Integer.parseInt(values[fieldsIndexMap.get("nextSpliceExonRank")]))
                    .nextSpliceExonPhase(Integer.parseInt(values[fieldsIndexMap.get("nextSpliceExonPhase")]))
                    .nextSpliceDistance(Integer.parseInt(values[fieldsIndexMap.get("nextSpliceDistance")]))
                    .totalExonCount(Integer.parseInt(values[fieldsIndexMap.get("totalExonCount")]))
                    .exonUp(Integer.parseInt(values[fieldsIndexMap.get("exonUp")]))
                    .exonDown(Integer.parseInt(values[fieldsIndexMap.get("exonDown")]))
                    .build());
        }

        return breakends;
    }

    private static String header()
    {
        return new StringJoiner(TSV_DELIM)
                .add("id")
                .add("svId")
                .add("vcfId")
                .add("coords")
                .add("isStart")
                .add("gene")
                .add("transcriptId")
                .add("canonical")
                .add("geneOrientation")
                .add("disruptive")
                .add("reportedStatus")
                .add("undisruptedCopyNumber")
                .add("regionType")
                .add("codingType")
                .add("biotype")
                .add("exonicBasePhase")
                .add("nextSpliceExonRank")
                .add("nextSpliceExonPhase")
                .add("nextSpliceDistance")
                .add("totalExonCount")
                .add("exonUp")
                .add("exonDown")
                .toString();
    }

    private static String toString(final LinxBreakend breakend)
    {
        return new StringJoiner(TSV_DELIM)
                .add(String.valueOf(breakend.id()))
                .add(String.valueOf(breakend.svId()))
                .add(breakend.vcfId())
                .add(breakend.coords())
                .add(String.valueOf(breakend.isStart()))
                .add(String.valueOf(breakend.gene()))
                .add(String.valueOf(breakend.transcriptId()))
                .add(String.valueOf(breakend.canonical()))
                .add(String.valueOf(breakend.geneOrientation()))
                .add(String.valueOf(breakend.disruptive()))
                .add(String.valueOf(breakend.reportedStatus()))
                .add(String.format("%.4f", breakend.undisruptedCopyNumber()))
                .add(String.valueOf(breakend.regionType()))
                .add(String.valueOf(breakend.codingType()))
                .add(String.valueOf(breakend.biotype()))
                .add(String.valueOf(breakend.exonicBasePhase()))
                .add(String.valueOf(breakend.nextSpliceExonRank()))
                .add(String.valueOf(breakend.nextSpliceExonPhase()))
                .add(String.valueOf(breakend.nextSpliceDistance()))
                .add(String.valueOf(breakend.totalExonCount()))
                .add(String.valueOf(breakend.exonUp()))
                .add(String.valueOf(breakend.exonDown()))
                .toString();
    }
}
