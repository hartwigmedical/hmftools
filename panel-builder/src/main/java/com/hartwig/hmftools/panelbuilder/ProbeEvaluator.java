package com.hartwig.hmftools.panelbuilder;

import static java.lang.Math.abs;
import static java.lang.String.format;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.panelbuilder.Utils.isDnaSequenceNormal;

import java.text.DecimalFormat;
import java.util.stream.Stream;

import org.jetbrains.annotations.NotNull;

// Common candidate probe filtering.
public class ProbeEvaluator
{
    public static Stream<Probe> evaluateProbes(Stream<Probe> probes, final Criteria criteria)
    {
        return probes.map(probe -> evaluateProbe(probe, criteria));
    }

    protected static Probe evaluateProbe(Probe probe, final Criteria criteria)
    {
        probe = probe.withEvalCriteria(criteria);

        probe = evaluateSequence(probe);
        if(probe.rejected())
        {
            return probe;
        }

        probe = evaluateGcContent(probe);
        if(probe.rejected())
        {
            return probe;
        }

        probe = evaluateQualityScore(probe);
        return probe;
    }

    private static Probe evaluateQualityScore(Probe probe)
    {
        Criteria criteria = requireNonNull(probe.evalCriteria());
        double qualityScore = requireNonNull(probe.qualityScore());
        if(!(qualityScore >= criteria.qualityScoreMin()))
        {
            probe = probe.withRejectionReason("QS");
        }
        return probe;
    }

    private static Probe evaluateSequence(Probe probe)
    {
        String sequence = requireNonNull(probe.sequence());
        if(!isDnaSequenceNormal(sequence))
        {
            probe = probe.withRejectionReason("sequence");
        }
        return probe;
    }

    private static Probe evaluateGcContent(Probe probe)
    {
        Criteria criteria = requireNonNull(probe.evalCriteria());
        double gcContent = requireNonNull(probe.gcContent());
        if(!(abs(gcContent - criteria.gcContentTarget()) <= criteria.gcContentTolerance()))
        {
            probe = probe.withRejectionReason("GC");
        }
        return probe;
    }

    public record Criteria(
            // Quality score must be >= this value.
            double qualityScoreMin,
            // Target GC content.
            double gcContentTarget,
            // How much +/- gcContentTarget to accept.
            double gcContentTolerance
    )
    {
        private static final DecimalFormat DECIMAL_FORMAT;

        static
        {
            DECIMAL_FORMAT = new DecimalFormat();
            DECIMAL_FORMAT.setMinimumFractionDigits(0);
        }

        public Criteria
        {
            if(!(qualityScoreMin > 0 && qualityScoreMin <= 1))
            {
                // Note quality score is always required, quality=0 is never acceptable.
                throw new IllegalArgumentException("qualityScoreMin must be in range (0, 1]");
            }
            if(!(gcContentMax() >= 0 && gcContentMin() <= 1))
            {
                throw new IllegalArgumentException("GC content range must overlap range [0, 1]");
            }
        }

        public double gcContentMin()
        {
            return gcContentTarget - gcContentTolerance;
        }

        public double gcContentMax()
        {
            return gcContentTarget + gcContentTolerance;
        }

        @NotNull
        @Override
        public String toString()
        {
            String str = format("QS>=%s", DECIMAL_FORMAT.format(qualityScoreMin));
            if(gcContentMin() > 0 || gcContentMax() < 1)
            {
                // Only show GC criteria if it does anything.
                str += format(" GC=%s+-%s", DECIMAL_FORMAT.format(gcContentTarget), DECIMAL_FORMAT.format(gcContentTolerance));
            }
            return str;
        }
    }
}
