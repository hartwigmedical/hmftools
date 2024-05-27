package com.hartwig.hmftools.sage.pipeline;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.SamSlicerFactory;
import com.hartwig.hmftools.sage.evidence.EvidenceStats;
import com.hartwig.hmftools.sage.evidence.ReadContextEvidence;
import com.hartwig.hmftools.sage.phase.VariantPhaser;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.bqr.BqrRecordMap;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.evidence.ReadContextCounters;
import com.hartwig.hmftools.sage.quality.MsiJitterCalcs;

public class EvidenceStage
{
    private final SageConfig mConfig;
    private final SamSlicerFactory mSamSlicerFactory;

    private final ReadContextEvidence mReadContextEvidence;
    private final VariantPhaser mVariantPhaser;

    public EvidenceStage(
            final SageConfig config, final RefGenomeInterface refGenome, final Map<String, BqrRecordMap> qualityRecalibrationMap,
            final MsiJitterCalcs msiJitterCalcs, final PhaseSetCounter phaseSetCounter, final SamSlicerFactory samSlicerFactory)
    {
        mConfig = config;
        mSamSlicerFactory = samSlicerFactory;

        mReadContextEvidence = new ReadContextEvidence(config, refGenome, qualityRecalibrationMap, msiJitterCalcs);
        mVariantPhaser = new VariantPhaser(phaseSetCounter);
    }

    public ReadContextCounters findEvidence(
            final ChrBaseRegion region, final String sampleType, final List<String> samples, final List<Candidate> candidates, boolean checkPhasing)
    {
        // search BAMs for evidence of each candidate variant
        if(samples.isEmpty())
            return new ReadContextCounters(mConfig, candidates);

        int sampleCount = samples.size();
        final ReadContextCounters readContextCounters = new ReadContextCounters(mConfig, candidates);

        for(int i = 0; i < samples.size(); i++)
        {
            final String sample = samples.get(i);

            boolean collectPhasingGroups = checkPhasing && (i == 0);

            List<ReadContextCounter> readCounters = mReadContextEvidence.collectEvidence(
                    candidates, sample, mSamSlicerFactory, collectPhasingGroups ? mVariantPhaser : null);

            readContextCounters.addCounters(readCounters, sampleCount);
        }

        SG_LOGGER.trace("region({}) gathered {} evidence for {} variants",
                region, sampleType, readContextCounters.candidateCount());

        return readContextCounters;
    }

    public VariantPhaser getVariantPhaser() { return mVariantPhaser; }

    public final int[] getSyncCounts() { return mReadContextEvidence.getSynCounts(); }
    public EvidenceStats getEvidenceStats() { return mReadContextEvidence.evidenceStats(); }
}
