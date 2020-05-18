package com.hartwig.hmftools.svtools.germline;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStrList;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.seIndex;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantImpl;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantLegImpl;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.common.fusion.GeneAnnotation;

import org.jetbrains.annotations.NotNull;

public class GermlineSV
{
    public final String SampleId;
    public final String Id;
    public final String Filter;
    public final String GridssFilter;
    public final double QualScore;
    public final StructuralVariantType Type;
    public final String ChrStart;
    public final String ChrEnd;
    public final long PosStart;
    public final long PosEnd;
    public final int OrientStart;
    public final int OrientEnd;
    public final int NormalREF;
    public final int NormalRP;
    public final double NormalRPQ;
    public final int NormalSR;
    public final double NormalSRQ;
    public final int NormalVF;
    public final String InsertSequence;
    public final String Homology;

    private StructuralVariant mSV;

    private final String mAffectedGenes;

    private String[] mAsmbSvIds;

    private final List<List<GeneAnnotation>> mBreakendGenes; // genes affected by start and end breakends
    private final List<GeneAnnotation> mOverlappedGenes; // genes fully overlapped by SV
    private final Map<GeneAnnotation,String> mDisruptions; // gene and the type of disruption

    private static final String DELIMITER = ",";
    public static final int FIELD_COUNT = 25;

    public GermlineSV(
            final String sampleId, final String id, final String filter, final String gridssFilter, double qualScore,
            final StructuralVariantType type, final String chrStart, final String chrEnd, long posStart, long posEnd,
            int orientStart, int orientEnd, int normalREF, int normalRP, double normalRPQ,
            int normalSR, double normalSRQ, int normalVF, final String insertSequence, final String homology,
            final String affectedGenes, final String asmbStart, final String asmbEnd)
    {
        SampleId = sampleId;
        Id = id;
        Filter = filter;
        GridssFilter = gridssFilter;
        QualScore = qualScore;
        Type = type;
        ChrStart = chrStart;
        ChrEnd = chrEnd;
        PosStart = posStart;
        PosEnd = posEnd;
        OrientStart = orientStart;
        OrientEnd = orientEnd;
        NormalREF = normalREF;
        NormalRP = normalRP;
        NormalRPQ = normalRPQ;
        NormalSR = normalSR;
        NormalSRQ = normalSRQ;
        NormalVF = normalVF;
        InsertSequence = insertSequence;
        Homology = homology;
        mAffectedGenes = affectedGenes;

        mAsmbSvIds = new String[SE_PAIR];
        mAsmbSvIds[SE_START] = asmbStart;
        mAsmbSvIds[SE_END] = asmbEnd;

        mSV = null;
        mBreakendGenes = Lists.newArrayListWithExpectedSize(2);
        mBreakendGenes.add(Lists.newArrayList());
        mBreakendGenes.add(Lists.newArrayList());

        mOverlappedGenes = Lists.newArrayList();

        mDisruptions = Maps.newHashMap();
    }

    public String affectedGenes()
    {
        String affectedGenes = "";

        final List<String> uniqueGenes = mOverlappedGenes.stream().map(x -> x.GeneName).collect(toList());

        for(int se = SE_START; se <= SE_END; ++se)
        {
            mBreakendGenes.get(se).stream().filter(x -> !uniqueGenes.contains(x.GeneName)).forEach(x -> uniqueGenes.add(x.GeneName));
        }

        return appendStrList(uniqueGenes, ';');
    }

    public String assemblySvIds(boolean isStart) { return mAsmbSvIds[seIndex(isStart)]; }
    public void setAssemblySvId(int seIndex, final String svId) { mAsmbSvIds[seIndex] = svId; }

    public StructuralVariant sv() { return mSV; }

    public final List<List<GeneAnnotation>> getBreakendGenes() { return mBreakendGenes; };
    public final List<GeneAnnotation> getOverlappedGenes() { return mOverlappedGenes; };
    public final List<GeneAnnotation> getBreakendGenes(int seIndex) { return mBreakendGenes.get(seIndex); };
    public final Map<GeneAnnotation,String> getDisruptions() { return mDisruptions; }

    private static final String FILE_EXTENSION = ".linx.germline_sv.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<GermlineSV> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<GermlineSV> svList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(svList));
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<GermlineSV> svList)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        svList.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static List<GermlineSV> fromLines(@NotNull List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith("SampleId")).map(x -> fromString(x, false)).collect(toList());
    }

    public static String stripBam(final String sampleId)
    {
        return sampleId.replaceAll("_dedup.realigned.bam","")
                .replaceAll(".sorted", "")
                .replaceAll(".bam", "");
    }

    @NotNull
    public static String header()
    {
        return new StringJoiner(DELIMITER)
            .add("SampleId")
            .add("Id")
            .add("Filter")
            .add("GridssFilter")
            .add("QualScore")
            .add("Type")
            .add("ChrStart")
            .add("ChrEnd")
            .add("PosStart")
            .add("PosEnd")
            .add("OrientStart")
            .add("OrientEnd")
            .add("NormalREF")
            .add("NormalRP")
            .add("NormalRPQ")
            .add("NormalSR")
            .add("NormalSRQ")
            .add("NormalVF")
            .add("InsertSequence")
            .add("Homology")
            .add("AffectedGenes")
            .add("AsmbStart")
            .add("AsmbEnd")
            .toString();
    }

    @NotNull
    public static String toString(@NotNull final GermlineSV sv)
    {
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(sv.SampleId))
                .add(String.valueOf(sv.Id))
                .add(String.valueOf(sv.Filter))
                .add(String.valueOf(sv.GridssFilter))
                .add(String.format("%.1f", sv.QualScore))
                .add(String.valueOf(sv.Type))
                .add(String.valueOf(sv.ChrStart))
                .add(String.valueOf(sv.ChrEnd))
                .add(String.valueOf(sv.PosStart))
                .add(String.valueOf(sv.PosEnd))
                .add(String.valueOf(sv.OrientStart))
                .add(String.valueOf(sv.OrientEnd))
                .add(String.valueOf(sv.NormalREF))
                .add(String.valueOf(sv.NormalRP))
                .add(String.format("%.2f", sv.NormalRPQ))
                .add(String.valueOf(sv.NormalSR))
                .add(String.format("%.2f", sv.NormalSRQ))
                .add(String.valueOf(sv.NormalVF))
                .add(String.valueOf(sv.InsertSequence))
                .add(String.valueOf(sv.Homology))
                .add(String.valueOf(sv.affectedGenes()))
                .add(String.valueOf(sv.assemblySvIds(true)))
                .add(String.valueOf(sv.assemblySvIds(false)))
                .toString();
    }

    @NotNull
    public static GermlineSV fromString(@NotNull final String sv, boolean reprocess)
    {
        String[] values = sv.split(DELIMITER, -1);
        int index = 0;

        return new GermlineSV(
                stripBam(values[index++]),
                values[index++],
                values[index++],
                values[index++],
                Double.parseDouble(values[index++]),
                StructuralVariantType.valueOf(values[index++]),
                values[index++],
                values[index++],
                Long.parseLong(values[index++]),
                Long.parseLong(values[index++]),
                Integer.parseInt(values[index++]),
                Integer.parseInt(values[index++]),
                Integer.parseInt(values[index++]),
                Integer.parseInt(values[index++]),
                Double.parseDouble(values[index++]),
                Integer.parseInt(values[index++]),
                Double.parseDouble(values[index++]),
                Integer.parseInt(values[index++]),
                values[index++],
                values[index++],
                reprocess ? "" : values[index++],
                reprocess ? "" : values[index++],
                reprocess ? "" : values[index++]);
    }

    public void createSV()
    {
        StructuralVariantLeg start = ImmutableStructuralVariantLegImpl.builder()
                .chromosome(ChrStart)
                .position(PosStart)
                .orientation((byte)OrientStart)
                .homology(Homology)
                .anchoringSupportDistance(0)
                .build();

        StructuralVariantLeg end = Type != SGL ? ImmutableStructuralVariantLegImpl.builder()
                .chromosome(ChrEnd)
                .position(PosEnd)
                .orientation((byte)OrientEnd)
                .homology(Homology)
                .anchoringSupportDistance(0)
                .build() : null;

        mSV = ImmutableStructuralVariantImpl.builder()
                .id(Id)
                .type(Type)
                .start(start)
                .end(end)
                .qualityScore(QualScore)
                .insertSequence(InsertSequence)
                .recovered(false)
                .build();
    }

    public String toString()
    {
        return String.format("%s:%s %s:%d - %s:%d", Id, Type.toString(), ChrStart, PosStart, ChrEnd, PosEnd);
    }
}
