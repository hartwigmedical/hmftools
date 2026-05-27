package com.hartwig.hmftools.panelbuilder;

import static java.lang.Math.abs;
import static java.lang.String.format;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.panelbuilder.SequenceUtils.buildSequence;

import java.text.DecimalFormat;
import java.util.function.Function;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import org.jetbrains.annotations.NotNull;

// Common candidate probe filtering.
public class ProbeEvaluator
{
    private final Function<Probe, SequenceData> mGetSequence;
    private final Function<Stream<Probe>, Stream<Probe>> mAnnotateQualityScore;

    protected ProbeEvaluator(final Function<Probe, SequenceData> getSequence,
            final Function<Stream<Probe>, Stream<Probe>> annotateQualityScore)
    {
        mGetSequence = getSequence;
        mAnnotateQualityScore = annotateQualityScore;
    }

    public ProbeEvaluator(final RefGenomeInterface refGenome, final ProbeQualityScorer probeQualityScorer)
    {
        this(
                probe -> buildSequence(refGenome, probe.definition()),
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
        probes = probes.map(this::annotateAndEvaluateSequenceAndGc);
        probes = mAnnotateQualityScore.apply(probes)
                .map(ProbeEvaluator::evaluateQualityScore)
                .map(ProbeEvaluator::acceptIfNotRejected);
        return probes;
    }

    private Probe annotateAndEvaluateSequenceAndGc(Probe probe)
    {
        SequenceData sequenceData = mGetSequence.apply(probe);

        String sequence = sequenceData.baseString();
        probe = probe.withSequence(sequence);
        if(!sequenceData.isNormal())
        {
            return probe.withEvaluationResult(EvaluationResult.reject("sequence"));
        }

        double gcContent = sequenceData.gcContent();
        probe = probe.withGcContent(gcContent);
        Criteria criteria = requireNonNull(probe.evaluationCriteria());
        if(!(abs(gcContent - criteria.gcContentTarget()) <= criteria.gcContentTolerance()))
        {
            return probe.withEvaluationResult(EvaluationResult.reject("GC"));
        }

        return probe;
    }

    private static Probe evaluateQualityScore(Probe probe)
    {
        if(!probe.rejected())
        {
            Criteria criteria = requireNonNull(probe.evaluationCriteria());
            double qualityScore = requireNonNull(probe.qualityScore());
            if(!(qualityScore >= criteria.qualityScoreMin()))
            {
                return probe.withEvaluationResult(EvaluationResult.reject("QS"));
            }
        }
        return probe;
    }

    private static Probe acceptIfNotRejected(Probe probe)
    {
        if(probe.rejected())
        {
            return probe;
        }
        else
        {
            return probe.withEvaluationResult(EvaluationResult.accept());
        }
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
