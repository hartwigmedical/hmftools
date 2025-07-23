package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;

import java.text.DecimalFormat;
import java.util.Objects;
import java.util.function.Consumer;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.refgenome.CachedRefGenome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

// TODO: unit test

// Common candidate probe evaluation and filtering.
public class ProbeEvaluator
{
    private final RefGenomeInterface mRefGenome;
    private final ProbeQualityProfile mQualityProfile;
    // Hook to catch all candidate probes for output.
    @Nullable
    private final Consumer<Probe> mCandidateCallback;

    private static final Logger LOGGER = LogManager.getLogger(ProbeEvaluator.class);

    public ProbeEvaluator(final RefGenomeInterface refGenome, final ProbeQualityProfile qualityProfile,
            final @Nullable Consumer<Probe> candidateCallback)
    {
        // During probe generation it's common to evaluate many nearby probes, which can be exploited with caching to improve performance.
        mRefGenome = refGenome instanceof CachedRefGenome ? refGenome : new CachedRefGenome(refGenome);
        mQualityProfile = qualityProfile;
        mCandidateCallback = candidateCallback;
    }

    public Probe evaluateProbe(final Probe probe, final Criteria criteria)
    {
        return evaluateProbe(probe.withEvalCriteria(criteria));
    }

    private Probe evaluateProbe(Probe probe)
    {
        probe = evaluateQualityScore(probe);
        if(!probe.rejected())
        {
            probe = evaluateGcContent(probe);
        }
        logCandidateProbe(probe);
        return probe;
    }

    public Stream<Probe> evaluateProbes(Stream<Probe> probes, final Criteria criteria)
    {
        return probes.map(probe -> evaluateProbe(probe, criteria));
    }

    private Probe evaluateQualityScore(Probe probe)
    {
        Criteria criteria = Objects.requireNonNull(probe.evalCriteria());
        double qualityScore = getProbeQuality(probe);
        probe = probe.withQualityScore(qualityScore);
        if(!(qualityScore >= criteria.qualityScoreMin()))
        {
            probe = probe.withRejectionReason("quality");
        }
        return probe;
    }

    private Probe evaluateGcContent(Probe probe)
    {
        Criteria criteria = Objects.requireNonNull(probe.evalCriteria());
        String sequence = getProbeSequence(probe);
        probe = probe.withSequence(sequence);
        double gcContent = calcGcPercent(sequence);
        probe = probe.withGcContent(gcContent);
        if(!(gcContent >= criteria.gcContentMin() && gcContent <= criteria.gcContentMax()))
        {
            probe = probe.withRejectionReason("gc");
        }
        return probe;
    }

    private double getProbeQuality(final Probe probe)
    {
        return mQualityProfile.computeQualityScore(probe.region()).orElseGet(() ->
        {
            // Never want to accept a probe with no quality score, so just return 0 in that case to simplify the code elsewhere.
            // Maybe be interesting to know when this happens because the probe quality profile ideally covers the whole genome.
            LOGGER.trace("Candidate probe not covered by probe quality profile so assuming qualityScore=0 probe={}", probe);
            return 0d;
        });
    }

    private String getProbeSequence(final Probe probe)
    {
        ChrBaseRegion region = probe.region();
        return mRefGenome.getBaseString(region.chromosome(), region.start(), region.end());
    }

    private void logCandidateProbe(final Probe probe)
    {
        LOGGER.trace("Evaluated probe: {}", probe);
        if(mCandidateCallback != null)
        {
            mCandidateCallback.accept(probe);
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
            String str = format("quality>=%s", DECIMAL_FORMAT.format(qualityScoreMin));
            if(gcContentMin() > 0 || gcContentMax() < 1)
            {
                // Only show GC criteria if it does anything.
                str += format(" gc=%s+-%s", DECIMAL_FORMAT.format(gcContentTarget), DECIMAL_FORMAT.format(gcContentTolerance));
            }
            return str;
        }
    }
}
