package com.hartwig.hmftools.svtools.germline;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.seIndex;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.ImmutableStructuralVariantImpl;
import com.hartwig.hmftools.common.sv.ImmutableStructuralVariantLegImpl;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

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

    private String[] mAsmbSvIds;

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

        mAsmbSvIds = new String[SE_PAIR];
        mAsmbSvIds[SE_START] = asmbStart;
        mAsmbSvIds[SE_END] = asmbEnd;

        mSV = null;
    }

    public String assemblySvIds(boolean isStart) { return mAsmbSvIds[seIndex(isStart)]; }
    public void setAssemblySvId(int seIndex, final String svId) { mAsmbSvIds[seIndex] = svId; }

    public StructuralVariant sv() { return mSV; }

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
                .hotspot(false)
                .build();
    }

    public String toString()
    {
        return String.format("%s:%s %s:%d - %s:%d", Id, Type.toString(), ChrStart, PosStart, ChrEnd, PosEnd);
    }
}
