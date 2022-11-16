package com.hartwig.hmftools.common.linx;

import static com.hartwig.hmftools.common.linx.LinxCluster.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.TranscriptCodingType;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.immutables.value.Value;

@Value.Immutable
public abstract class LinxBreakend
{
    public abstract int id();
    public abstract int svId();
    public abstract boolean isStart();
    public abstract String gene();
    public abstract String transcriptId();
    public abstract boolean canonical();
    public abstract String geneOrientation();
    public abstract boolean disruptive();
    public abstract boolean reportedDisruption();
    public abstract double undisruptedCopyNumber();
    public abstract TranscriptRegionType regionType();
    public abstract TranscriptCodingType codingType();
    public abstract String biotype();
    public abstract int exonicBasePhase();
    public abstract int nextSpliceExonRank();
    public abstract int nextSpliceExonPhase();
    public abstract int nextSpliceDistance();
    public abstract int totalExonCount();

    // additional fields for patient report
    public abstract StructuralVariantType type();
    public abstract String chromosome();
    public abstract int orientation();
    public abstract int strand();
    public abstract String chrBand();
    public abstract int exonUp();
    public abstract int exonDown();
    public abstract double junctionCopyNumber();

    private static final String FILE_EXTENSION = ".linx.breakend.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static List<LinxBreakend> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(final String filename, List<LinxBreakend> breakends) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(breakends));
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
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIMITER);
        lines.remove(0);

        List<LinxBreakend> breakends = Lists.newArrayList();

        // backwards compatibility on old column name
        Integer codingTypeIndex = fieldsIndexMap.containsValue("codingType") ?
                fieldsIndexMap.get("codingType") : fieldsIndexMap.get("codingContext");

        for(int i = 0; i < lines.size(); ++i)
        {
            String[] values = lines.get(i).split(DELIMITER);

            breakends.add(ImmutableLinxBreakend.builder()
                    .id(Integer.parseInt(values[fieldsIndexMap.get("id")]))
                    .svId(Integer.parseInt(values[fieldsIndexMap.get("svId")]))
                    .isStart(Boolean.parseBoolean(values[fieldsIndexMap.get("isStart")]))
                    .gene(values[fieldsIndexMap.get("gene")])
                    .transcriptId(values[fieldsIndexMap.get("transcriptId")])
                    .canonical(Boolean.parseBoolean(values[fieldsIndexMap.get("canonical")]))
                    .geneOrientation(values[fieldsIndexMap.get("geneOrientation")])
                    .disruptive(Boolean.parseBoolean(values[fieldsIndexMap.get("disruptive")]))
                    .reportedDisruption(Boolean.parseBoolean(values[fieldsIndexMap.get("reportedDisruption")]))
                    .undisruptedCopyNumber(Double.parseDouble(values[fieldsIndexMap.get("undisruptedCopyNumber")]))
                    .regionType(TranscriptRegionType.valueOf(values[fieldsIndexMap.get("regionType")]))
                    .codingType(TranscriptCodingType.valueOf(values[codingTypeIndex]))
                    .biotype(values[fieldsIndexMap.get("biotype")])
                    .exonicBasePhase(Integer.parseInt(values[fieldsIndexMap.get("exonicBasePhase")]))
                    .nextSpliceExonRank(Integer.parseInt(values[fieldsIndexMap.get("nextSpliceExonRank")]))
                    .nextSpliceExonPhase(Integer.parseInt(values[fieldsIndexMap.get("nextSpliceExonPhase")]))
                    .nextSpliceDistance(Integer.parseInt(values[fieldsIndexMap.get("nextSpliceDistance")]))
                    .totalExonCount(Integer.parseInt(values[fieldsIndexMap.get("totalExonCount")]))
                    .type(StructuralVariantType.valueOf(values[fieldsIndexMap.get("type")]))
                    .chromosome(values[fieldsIndexMap.get("chromosome")])
                    .orientation(Integer.parseInt(values[fieldsIndexMap.get("orientation")]))
                    .strand(Integer.parseInt(values[fieldsIndexMap.get("strand")]))
                    .chrBand(values[fieldsIndexMap.get("chrBand")])
                    .exonUp(Integer.parseInt(values[fieldsIndexMap.get("exonUp")]))
                    .exonDown(Integer.parseInt(values[fieldsIndexMap.get("exonDown")]))
                    .junctionCopyNumber(Double.parseDouble(values[fieldsIndexMap.get("junctionCopyNumber")]))
                    .build());
        }

        return breakends;
    }

    private static String header()
    {
        return new StringJoiner(DELIMITER)
                .add("id")
                .add("svId")
                .add("isStart")
                .add("gene")
                .add("transcriptId")
                .add("canonical")
                .add("geneOrientation")
                .add("disruptive")
                .add("reportedDisruption")
                .add("undisruptedCopyNumber")
                .add("regionType")
                .add("codingType")
                .add("biotype")
                .add("exonicBasePhase")
                .add("nextSpliceExonRank")
                .add("nextSpliceExonPhase")
                .add("nextSpliceDistance")
                .add("totalExonCount")
                .add("type")
                .add("chromosome")
                .add("orientation")
                .add("strand")
                .add("chrBand")
                .add("exonUp")
                .add("exonDown")
                .add("junctionCopyNumber")
                .toString();
    }

    private static String toString(final LinxBreakend breakend)
    {
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(breakend.id()))
                .add(String.valueOf(breakend.svId()))
                .add(String.valueOf(breakend.isStart()))
                .add(String.valueOf(breakend.gene()))
                .add(String.valueOf(breakend.transcriptId()))
                .add(String.valueOf(breakend.canonical()))
                .add(String.valueOf(breakend.geneOrientation()))
                .add(String.valueOf(breakend.disruptive()))
                .add(String.valueOf(breakend.reportedDisruption()))
                .add(String.format("%.4f", breakend.undisruptedCopyNumber()))
                .add(String.valueOf(breakend.regionType()))
                .add(String.valueOf(breakend.codingType()))
                .add(String.valueOf(breakend.biotype()))
                .add(String.valueOf(breakend.exonicBasePhase()))
                .add(String.valueOf(breakend.nextSpliceExonRank()))
                .add(String.valueOf(breakend.nextSpliceExonPhase()))
                .add(String.valueOf(breakend.nextSpliceDistance()))
                .add(String.valueOf(breakend.totalExonCount()))
                .add(String.valueOf(breakend.type()))
                .add(String.valueOf(breakend.chromosome()))
                .add(String.valueOf(breakend.orientation()))
                .add(String.valueOf(breakend.strand()))
                .add(String.valueOf(breakend.chrBand()))
                .add(String.valueOf(breakend.exonUp()))
                .add(String.valueOf(breakend.exonDown()))
                .add(String.format("%.4f", breakend.junctionCopyNumber()))
                .toString();
    }
}
