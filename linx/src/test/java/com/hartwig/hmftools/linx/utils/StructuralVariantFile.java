package com.hartwig.hmftools.linx.utils;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.sv.SvVcfTags.INFERRED;

import com.hartwig.hmftools.common.sv.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

public class StructuralVariantFile
{
    public static StructuralVariantData fromString(final String[] values)
    {
        int index = 0;

        final ImmutableStructuralVariantData.Builder builder = ImmutableStructuralVariantData.builder()
                .id(Integer.parseInt(values[index++]))
                .vcfIdStart(values[index++])
                .vcfIdEnd(values[index++])
                .startChromosome(values[index++])
                .endChromosome(getStrValue(values[index++]))
                .startPosition(Integer.parseInt(values[index++]))
                .endPosition(getIntValue(values[index++]))
                .startOrientation(Byte.parseByte(values[index++]))
                .endOrientation(getByteValue(values[index++]))
                .startHomologySequence(values[index++])
                .endHomologySequence(getStrValue(values[index++]))
                .startAF(Double.parseDouble(values[index++]))
                .endAF(getDoubleValue(values[index++]))
                .junctionCopyNumber(Double.parseDouble(values[index++]))
                .adjustedStartAF(Double.parseDouble(values[index++]))
                .adjustedEndAF(getDoubleValue(values[index++]))
                .adjustedStartCopyNumber(Double.parseDouble(values[index++]))
                .adjustedEndCopyNumber(getDoubleValue(values[index++]))
                .adjustedStartCopyNumberChange(Double.parseDouble(values[index++]))
                .adjustedEndCopyNumberChange(getDoubleValue(values[index++]))
                .insertSequence(values[index++]);

        StructuralVariantType type = StructuralVariantType.valueOf(values[index++]);
        final String filterStr = values[index++];
        if(type == SGL && filterStr.equals(INFERRED))
            type = INF;

        builder.type(type)
                .filter(filterStr)
                .qualityScore(Double.parseDouble(values[index++]))
                .event(values[index++])
                .startTumorVariantFragmentCount(Integer.parseInt(values[index++]))
                .startTumorReferenceFragmentCount(Integer.parseInt(values[index++]))
                .startNormalVariantFragmentCount(Integer.parseInt(values[index++]))
                .startNormalReferenceFragmentCount(Integer.parseInt(values[index++]))
                .endTumorVariantFragmentCount(getIntValue(values[index++]))
                .endTumorReferenceFragmentCount(getIntValue(values[index++]))
                .endNormalVariantFragmentCount(getIntValue(values[index++]))
                .endNormalReferenceFragmentCount(getIntValue(values[index++]))
                .startIntervalOffsetStart(Integer.parseInt(values[index++]))
                .startIntervalOffsetEnd(Integer.parseInt(values[index++]))
                .endIntervalOffsetStart(getIntValue(values[index++]))
                .endIntervalOffsetEnd(getIntValue(values[index++]))
                .inexactHomologyOffsetStart(Integer.parseInt(values[index++]))
                .inexactHomologyOffsetEnd(getIntValue(values[index++]))
                .startLinkedBy(values[index++])
                .endLinkedBy(values[index++])
                .insertSequenceAlignments(values[index++])
                .insertSequenceRepeatClass(values[index++])
                .insertSequenceRepeatType(values[index++])
                .insertSequenceRepeatOrientation(getByteValue(values[index++]))
                .insertSequenceRepeatCoverage(getDoubleValue(values[index++]))
                .startAnchoringSupportDistance(Integer.parseInt(values[index++]))
                .endAnchoringSupportDistance(getIntValue(values[index++]))
                .ponCount(0);

        return builder.build();
    }

    private static String getStrValue(final String value) { return value.equals("NULL") ? "" : value; }
    private static double getDoubleValue(final String value) { return value.equals("NULL") ? 0 : Double.parseDouble(value); }
    private static int getIntValue(final String value) { return value.equals("NULL") ? 0 : Integer.parseInt(value); }
    private static long getLongValue(final String value) { return value.equals("NULL") ? 0 : Long.parseLong(value); }
    private static byte getByteValue(final String value) { return value.equals("NULL") ? 0 : Byte.valueOf(value); }
}
