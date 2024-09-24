package com.hartwig.hmftools.common.linx;

import static java.lang.String.valueOf;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getDoubleValue;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getIntValue;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getValue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

public final class LinxGermlineDisruption
{
    public final int SvId;
    public final String VcfId;
    public final String ChromosomeStart;
    public final String ChromosomeEnd;
    public final int PositionStart;
    public final int PositionEnd;
    public final byte OrientStart;
    public final byte OrientEnd;
    public final StructuralVariantType Type;
    public final String Filter;
    public final String EventId;
    public final double QualScore;
    public final String HomologyStart;
    public final String HomologyEnd;
    public final int InexactHomologyStart;
    public final int InexactHomologyEnd;
    public final double JunctionCopyNumber;
    public final double AdjustedAFStart;
    public final double AdjustedAFEnd;
    public final double AdjustedCopyNumberStart;
    public final double AdjustedCopyNumberEnd;
    public final double AdjustedCopyNumberChangeStart;
    public final double AdjustedCopyNumberChangeEnd;
    public final int GermlineFragments;
    public final int GermlineReferenceFragmentsStart;
    public final int GermlineReferenceFragmentsEnd;
    public final int TumorFragments;
    public final int TumorReferenceFragmentsStart;
    public final int TumorReferenceFragmentsEnd;
    public final String InsertSequence;
    public final String InsertSequenceAlignments;
    public final String InsertSequenceRepeatClass;
    public final String InsertSequenceRepeatType;
    public final String GeneName;
    public final int ClusterId;
    public final int ClusterCount;
    public final String ResolvedType;
    public final String LinkedByStart;
    public final String LinkedByEnd;
    public final int CohortFrequency;

    public LinxGermlineDisruption(
            final int svId, final String vcfId, final String chromosomeStart, final String chromosomeEnd,
            final int positionStart, final int positionEnd, final byte orientStart, final byte orientEnd, final StructuralVariantType type,
            final String filter, final String eventId, final double qualScore,
            final String homologyStart, final String homologyEnd, final int inexactHomologyStart, final int inexactHomologyEnd,
            final double junctionCopyNumber,
            final double adjAFStart, final double adjAFEnd, final double adjCopyNumberStart, final double adjCopyNumberEnd,
            final double adjustedCopyNumberChangeStart, final double adjustedCopyNumberChangeEnd,
            final int germlineFragments, final int germlineReferenceFragmentsStart, final int germlineReferenceFragmentsEnd,
            final int tumorFragments, final int tumorReferenceFragmentsStart, final int tumorReferenceFragmentsEnd,
            final String insSeq, final String insSeqAlignments, final String insSeqRepeatClass, final String insSeqRepeatType,
            final String geneName, final int clusterId, final int clusterCount, final String resolvedType,
            final String linkedByStart, final String linkedByEnd, final int cohortFrequency)
    {
        SvId = svId;
        VcfId = vcfId;
        ChromosomeStart = chromosomeStart;
        ChromosomeEnd = chromosomeEnd;
        PositionStart = positionStart;
        PositionEnd = positionEnd;
        OrientStart = orientStart;
        OrientEnd = orientEnd;
        Type = type;
        Filter = filter;
        EventId = eventId;
        QualScore = qualScore;
        HomologyStart = homologyStart;
        HomologyEnd = homologyEnd;
        InexactHomologyStart = inexactHomologyStart;
        InexactHomologyEnd = inexactHomologyEnd;
        JunctionCopyNumber = junctionCopyNumber;
        AdjustedAFStart = adjAFStart;
        AdjustedAFEnd = adjAFEnd;
        AdjustedCopyNumberStart = adjCopyNumberStart;
        AdjustedCopyNumberEnd = adjCopyNumberEnd;
        AdjustedCopyNumberChangeStart = adjustedCopyNumberChangeStart;
        AdjustedCopyNumberChangeEnd = adjustedCopyNumberChangeEnd;
        GermlineFragments = germlineFragments;
        GermlineReferenceFragmentsStart = germlineReferenceFragmentsStart;
        GermlineReferenceFragmentsEnd = germlineReferenceFragmentsEnd;
        TumorFragments = tumorFragments;
        TumorReferenceFragmentsStart = tumorReferenceFragmentsStart;
        TumorReferenceFragmentsEnd = tumorReferenceFragmentsEnd;
        InsertSequence = insSeq;
        InsertSequenceAlignments = insSeqAlignments;
        InsertSequenceRepeatClass = insSeqRepeatClass;
        InsertSequenceRepeatType = insSeqRepeatType;
        GeneName = geneName;
        ClusterId = clusterId;
        ClusterCount = clusterCount;
        ResolvedType = resolvedType;
        LinkedByStart = linkedByStart;
        LinkedByEnd = linkedByEnd;
        CohortFrequency = cohortFrequency;
    }

    private static final String EXTENSION = ".linx.germline.disruption.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static List<LinxGermlineDisruption> read(final String fileName) throws IOException
    {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    public static void write(final String fileName, List<LinxGermlineDisruption> deletions) throws IOException
    {
        Files.write(new File(fileName).toPath(), toLines(deletions));
    }

    private static List<String> toLines(final List<LinxGermlineDisruption> deletions)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        deletions.forEach(x -> lines.add(toString(x)));
        return lines;
    }

    private static String header()
    {
        return new StringJoiner(TSV_DELIM, "", "")
                .add("svId")
                .add("vcfId")
                .add("chromosomeStart")
                .add("chromosomeEnd")
                .add("positionStart")
                .add("positionEnd")
                .add("orientStart")
                .add("orientEnd")
                .add("type")
                .add("filter")
                .add("event")
                .add("qualScore")
                .add("homologySequenceStart")
                .add("homologySequenceEnd")
                .add("inexactHomologyStart")
                .add("inexactHomologyEnd")
                .add("junctionCopyNumber")
                .add("adjustedAFStart")
                .add("adjustedAFEnd")
                .add("adjustedCopyNumberStart")
                .add("adjustedCopyNumberEnd")
                .add("adjustedCopyNumberChangeStart")
                .add("adjustedCopyNumberChangeEnd")
                .add("germlineFragments")
                .add("germlineReferenceFragmentsStart")
                .add("germlineReferenceFragmentsEnd")
                .add("tumorFragments")
                .add("tumorReferenceFragmentsStart")
                .add("tumorReferenceFragmentsEnd")
                .add("insertSequence")
                .add("insertSequenceAlignments")
                .add("insertSequenceRepeatClass")
                .add("insertSequenceRepeatType")
                .add("gene")
                .add("clusterId")
                .add("clusterCount")
                .add("resolvedType")
                .add("linkedByStart")
                .add("linkedByEnd")
                .add("cohortFrequency")
                .toString();
    }

    private static String toString(final LinxGermlineDisruption germlineSv)
    {
        return new StringJoiner(TSV_DELIM)
                .add(valueOf(germlineSv.SvId))
                .add(germlineSv.VcfId)
                .add(germlineSv.ChromosomeStart)
                .add(germlineSv.ChromosomeEnd)
                .add(valueOf(germlineSv.PositionStart))
                .add(valueOf(germlineSv.PositionEnd))
                .add(valueOf(germlineSv.OrientStart))
                .add(valueOf(germlineSv.OrientEnd))
                .add(valueOf(germlineSv.Type))
                .add(germlineSv.Filter)
                .add(germlineSv.EventId)
                .add(valueOf(germlineSv.QualScore))
                .add(germlineSv.HomologyStart)
                .add(germlineSv.HomologyEnd)
                .add(valueOf(germlineSv.InexactHomologyStart))
                .add(valueOf(germlineSv.InexactHomologyEnd))
                .add(valueOf(germlineSv.JunctionCopyNumber))
                .add(valueOf(germlineSv.AdjustedAFStart))
                .add(valueOf(germlineSv.AdjustedAFEnd))
                .add(valueOf(germlineSv.AdjustedCopyNumberStart))
                .add(valueOf(germlineSv.AdjustedCopyNumberEnd))
                .add(valueOf(germlineSv.AdjustedCopyNumberChangeStart))
                .add(valueOf(germlineSv.AdjustedCopyNumberChangeEnd))
                .add(valueOf(germlineSv.GermlineFragments))
                .add(valueOf(germlineSv.GermlineReferenceFragmentsStart))
                .add(valueOf(germlineSv.GermlineReferenceFragmentsEnd))
                .add(valueOf(germlineSv.TumorFragments))
                .add(valueOf(germlineSv.TumorReferenceFragmentsStart))
                .add(valueOf(germlineSv.TumorReferenceFragmentsEnd))
                .add(germlineSv.InsertSequence)
                .add(germlineSv.InsertSequenceAlignments)
                .add(germlineSv.InsertSequenceRepeatClass)
                .add(germlineSv.InsertSequenceRepeatType)
                .add(germlineSv.GeneName)
                .add(valueOf(germlineSv.ClusterId))
                .add(valueOf(germlineSv.ClusterCount))
                .add(germlineSv.ResolvedType)
                .add(germlineSv.LinkedByStart)
                .add(germlineSv.LinkedByEnd)
                .add(valueOf(germlineSv.CohortFrequency))
                .toString();
    }

    static List<LinxGermlineDisruption> fromLines(final List<String> lines)
    {
        final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);
        lines.remove(0);

        List<LinxGermlineDisruption> germlineSVs = Lists.newArrayList();

        for(final String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            germlineSVs.add(new LinxGermlineDisruption(
                    getIntValue(fieldsIndexMap, "svId", values),
                    getValue(fieldsIndexMap, "vcfId", "", values),
                    values[fieldsIndexMap.get("chromosomeStart")],
                    values[fieldsIndexMap.get("chromosomeEnd")],
                    getIntValue(fieldsIndexMap, "positionStart", values),
                    getIntValue(fieldsIndexMap, "positionEnd", values),
                    Byte.parseByte(values[fieldsIndexMap.get("orientStart")]),
                    Byte.parseByte(values[fieldsIndexMap.get("orientEnd")]),
                    StructuralVariantType.valueOf(values[fieldsIndexMap.get("type")]),
                    values[fieldsIndexMap.get("filter")], values[fieldsIndexMap.get("event")],
                    getDoubleValue(fieldsIndexMap, "qualScore", values),
                    getValue(fieldsIndexMap, "homologySequenceStart", "", values),
                    getValue(fieldsIndexMap, "homologySequenceEnd", "", values),
                    getIntValue(fieldsIndexMap, "inexactHomologyStart", values),
                    getIntValue(fieldsIndexMap, "inexactHomologyEnd", values),
                    getDoubleValue(fieldsIndexMap, "junctionCopyNumber", 0, values),
                    getDoubleValue(fieldsIndexMap, "adjustedAFStart", 0, values),
                    getDoubleValue(fieldsIndexMap, "adjustedAFEnd", 0, values),
                    getDoubleValue(fieldsIndexMap, "adjustedCopyNumberStart", 0, values),
                    getDoubleValue(fieldsIndexMap, "adjustedCopyNumberEnd", 0, values),
                    getDoubleValue(fieldsIndexMap, "adjustedCopyNumberChangeStart", 0, values),
                    getDoubleValue(fieldsIndexMap, "adjustedCopyNumberChangeEnd", 0, values),
                    getIntValue(fieldsIndexMap, "germlineFragments", values),
                    getIntValue(fieldsIndexMap, "germlineReferenceFragmentsStart", values),
                    getIntValue(fieldsIndexMap, "germlineReferenceFragmentsEnd", values),
                    getIntValue(fieldsIndexMap, "tumorFragments", values),
                    getIntValue(fieldsIndexMap, "tumorReferenceFragmentsStart", values),
                    getIntValue(fieldsIndexMap, "tumorReferenceFragmentsEnd", values),
                    values[fieldsIndexMap.get("insertSequence")],
                    values[fieldsIndexMap.get("insertSequenceAlignments")],
                    values[fieldsIndexMap.get("insertSequenceRepeatClass")],
                    values[fieldsIndexMap.get("insertSequenceRepeatType")],
                    getValue(fieldsIndexMap, "gene", "", values),
                    getIntValue(fieldsIndexMap, "clusterId", values),
                    getIntValue(fieldsIndexMap, "clusterCount", values),
                    values[fieldsIndexMap.get("resolvedType")],
                    values[fieldsIndexMap.get("linkedByStart")],
                    values[fieldsIndexMap.get("linkedByEnd")],
                    getIntValue(fieldsIndexMap, "cohortFrequency", values)));
        }

        return germlineSVs;
    }
}
