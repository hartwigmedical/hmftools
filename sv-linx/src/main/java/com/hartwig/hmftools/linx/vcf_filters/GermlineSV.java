package com.hartwig.hmftools.linx.vcf_filters;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantImpl;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantLegImpl;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

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
    public final String GenesStart;
    public final String GenesEnd;
    public final String GenePanelOverlaps;
    public final String AsmbStart;
    public final String AsmbEnd;

    private boolean mDisruptive;

    private static final String DELIMITER = ",";
    public static final int FIELD_COUNT = 25;

    public GermlineSV(
            final String sampleId, final String id, final String filter, final String gridssFilter, double qualScore,
            final StructuralVariantType type, final String chrStart, final String chrEnd, long posStart, long posEnd,
            int orientStart, int orientEnd, int normalREF, int normalRP, double normalRPQ,
            int normalSR, double normalSRQ, int normalVF, final String insertSequence, final String homology,
            final String genesStart, final String genesEnd, final String genePanelOverlaps, final String asmbStart, final String asmbEnd)
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
        GenesStart = genesStart;
        GenesEnd = genesEnd;
        GenePanelOverlaps = genePanelOverlaps;
        AsmbStart = asmbStart;
        AsmbEnd = asmbEnd;

        mDisruptive = false;
    }

    public boolean disruptive() { return mDisruptive; }
    public void setDisruptive(boolean toggle) { mDisruptive = toggle; }

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
        return lines.stream().filter(x -> !x.startsWith("SampleId")).map(GermlineSV::fromString).collect(toList());
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
            .add("GenesStart")
            .add("GenesEnd")
            .add("GenePanelOverlaps")
            .add("AsmbStart")
            .add("AsmbEnd")
            .toString();
    }

    @NotNull
    private static String toString(@NotNull final GermlineSV sv)
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
                .add(String.valueOf(sv.GenesStart))
                .add(String.valueOf(sv.GenesEnd))
                .add(String.valueOf(sv.GenePanelOverlaps))
                .add(String.valueOf(sv.AsmbStart))
                .add(String.valueOf(sv.AsmbEnd))
                .toString();
    }

    @NotNull
    public static GermlineSV fromString(@NotNull final String sv)
    {
        String[] values = sv.split(DELIMITER, -1);
        int index = 0;

        return new GermlineSV(
                values[index++],
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
                values[index++],
                values[index++],
                values[index++],
                values[index++],
                values[index++]);
    }

    public static StructuralVariant convert(final GermlineSV sv)
    {
        StructuralVariantLeg start = ImmutableStructuralVariantLegImpl.builder()
                .chromosome(sv.ChrStart)
                .position(sv.PosStart)
                .orientation((byte)sv.OrientStart)
                .homology(sv.Homology)
                .anchoringSupportDistance(0)
                .build();

        StructuralVariantLeg end = sv.Type != SGL ? ImmutableStructuralVariantLegImpl.builder()
                .chromosome(sv.ChrEnd)
                .position(sv.PosEnd)
                .orientation((byte)sv.OrientEnd)
                .homology(sv.Homology)
                .anchoringSupportDistance(0)
                .build() : null;

        return ImmutableStructuralVariantImpl.builder()
                .id(sv.Id)
                .type(sv.Type)
                .start(start)
                .end(end)
                .qualityScore(sv.QualScore)
                .insertSequence(sv.InsertSequence)
                .recovered(false)
                .build();

    }

    public String toString()
    {
        return String.format("%s:%s %s:%d - %s:%d", Id, Type.toString(), ChrStart, PosStart, ChrEnd, PosEnd);
    }
}
