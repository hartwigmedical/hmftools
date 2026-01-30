package com.hartwig.hmftools.panelbuilder;

import static java.lang.Math.abs;
import static java.lang.String.format;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;
import static com.hartwig.hmftools.panelbuilder.SequenceUtils.buildSequence;
import static com.hartwig.hmftools.panelbuilder.SequenceUtils.isDnaSequenceNormal;

import java.text.DecimalFormat;
import java.util.function.Function;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import org.jetbrains.annotations.NotNull;

// Common candidate probe filtering.
public class ProbeEvaluator
{
    private final Function<Probe, Probe> mAnnotateSequence;
    private final Function<Probe, Probe> mAnnotateGcContent;
    private final Function<Stream<Probe>, Stream<Probe>> mAnnotateQualityScore;

    protected ProbeEvaluator(final Function<Probe, Probe> annotateSequence, final Function<Probe, Probe> annotateGcContent,
            final Function<Stream<Probe>, Stream<Probe>> annotateQualityScore)
    {
        mAnnotateSequence = annotateSequence;
        mAnnotateGcContent = annotateGcContent;
        mAnnotateQualityScore = annotateQualityScore;
    }

    public ProbeEvaluator(final RefGenomeInterface refGenome, final ProbeQualityScorer probeQualityScorer)
    {
        this(
                probe -> probe.withSequence(buildSequence(refGenome, probe.definition())),
                probe -> probe.withGcContent(calcGcPercent(requireNonNull(probe.sequence()))),
                probeQualityScorer::computeQualityScores
        );
    }

    // Probes must have evaluation criteria already set.
    public Stream<Probe> evaluateProbes(Stream<Probe> probes)
    {
        // Efficient evaluation methodology:
        // 1. Annotate sequence
        // 2. Evaluate sequence
        // 3. If rejected, done
        // 4. Annotate GC
        // 5. Evaluate GC
        // 6. If rejected, done
        // 7. Try to annotate QS with probe quality profile
        // 8. If no QS, annotate QS with probe quality model
        // 9. Evaluate QS
        // 10. Done

        probes = probes.peek(probe ->
        {
            if(probe.evaluated())
            {
                throw new IllegalArgumentException("Probe must not already be evaluated");
            }
        });
        probes = probes
                .map(applyIfNotRejected(probe -> evaluateSequence(mAnnotateSequence.apply(probe))))
                // TODO? don't have to annotate GC if there's no limit
                .map(applyIfNotRejected(probe -> evaluateGcContent(mAnnotateGcContent.apply(probe))));
                // TODO? don't have to annotate QS if threshold is 0
        probes = mAnnotateQualityScore.apply(probes)
                .map(applyIfNotRejected(ProbeEvaluator::evaluateQualityScore));
        probes = probes.map(applyIfNotRejected(probe -> probe.withEvaluationResult(EvaluationResult.accept())));
        return probes;
    }

    private static Function<Probe, Probe> applyIfNotRejected(final Function<Probe, Probe> function)
    {
        return probe -> probe.rejected() ? probe : function.apply(probe);
    }

    private static Probe evaluateSequence(Probe probe)
    {
        String sequence = requireNonNull(probe.sequence());
        if(!isDnaSequenceNormal(sequence))
        {
            probe = probe.withEvaluationResult(EvaluationResult.reject("sequence"));
        }
        return probe;
    }

    private static Probe evaluateGcContent(Probe probe)
    {
        Criteria criteria = requireNonNull(probe.evaluationCriteria());
        double gcContent = requireNonNull(probe.gcContent());
        if(!(abs(gcContent - criteria.gcContentTarget()) <= criteria.gcContentTolerance()))
        {
            probe = probe.withEvaluationResult(EvaluationResult.reject("GC"));
        }
        return probe;
    }

    private static Probe evaluateQualityScore(Probe probe)
    {
        Criteria criteria = requireNonNull(probe.evaluationCriteria());
        double qualityScore = requireNonNull(probe.qualityScore());
        if(!(qualityScore >= criteria.qualityScoreMin()))
        {
            probe = probe.withEvaluationResult(EvaluationResult.reject("QS"));
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
            if(!(qualityScoreMin >= 0 && qualityScoreMin <= 1))
            {
                throw new IllegalArgumentException("qualityScoreMin must be in range [0, 1]");
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
