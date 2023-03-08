package com.hartwig.hmftools.common.aligner;

import java.util.Collection;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public enum AlignmentOperator
{
    MATCH('M'),
    MISMATCH('S'),
    INSERTION('+'),
    DELETION('-');

    public final char code;

    AlignmentOperator(char code)
    {
        this.code = code;
    }

    // convert the list of align ops to string like
    // MMMMSIDDMM
    @NotNull
    public static String toString(@NotNull Collection<AlignmentOperator> alignOps)
    {
        StringBuilder b = new StringBuilder();
        alignOps.stream().map(x -> x.code).forEach(b::append);
        return b.toString();
    }

    public static void logAlignment(
            @NotNull Logger logger,
            @NotNull Level logLevel,
            @NotNull String leftSeq, @NotNull String rightSeq,
            @NotNull Collection<AlignmentOperator> alignOps)
    {
        logAlignment(logger, logLevel, leftSeq, rightSeq, 0, 0, alignOps);
    }

    public static void logAlignment(
            @NotNull Logger logger,
            @NotNull Level logLevel,
            @NotNull String leftSeq, @NotNull String rightSeq,
            int leftSeqAlignStart, int rightSeqAlignStart,
            @NotNull Collection<AlignmentOperator> alignOps)
    {
        if (!logger.isEnabled(logLevel))
            return;

        StringBuilder seqAlignBuilder = new StringBuilder();
        StringBuilder alignOpBuilder = new StringBuilder();
        StringBuilder refSeqAlignBuilder = new StringBuilder();

        int i = leftSeqAlignStart;
        int j = rightSeqAlignStart;

        for (AlignmentOperator op : alignOps)
        {
            switch (op)
            {
                case MATCH:
                    seqAlignBuilder.append(leftSeq.charAt(i++));
                    refSeqAlignBuilder.append(rightSeq.charAt(j++));
                    alignOpBuilder.append('|');
                    break;
                case MISMATCH:
                    seqAlignBuilder.append(leftSeq.charAt(i++));
                    refSeqAlignBuilder.append(rightSeq.charAt(j++));
                    alignOpBuilder.append(' ');
                    break;
                case INSERTION:
                    seqAlignBuilder.append(leftSeq.charAt(i++));
                    refSeqAlignBuilder.append('-');
                    alignOpBuilder.append(' ');
                    break;
                case DELETION:
                    seqAlignBuilder.append('-');
                    refSeqAlignBuilder.append(rightSeq.charAt(j++));
                    alignOpBuilder.append(' ');
                    break;
            }
        }

        logger.log(logLevel, "alignment:");
        logger.log(logLevel, "{}", AlignmentOperator.toString(alignOps));
        logger.log(logLevel, "{}", seqAlignBuilder);
        logger.log(logLevel, "{}", alignOpBuilder);
        logger.log(logLevel, "{}", refSeqAlignBuilder);
    }
}
