package com.hartwig.hmftools.neo.score;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.FlankCounts.FLANK_AA_COUNT;
import static com.hartwig.hmftools.neo.score.DataLoader.loadAlleleCoverage;
import static com.hartwig.hmftools.neo.score.DataLoader.loadNeoEpitopes;
import static com.hartwig.hmftools.neo.score.DataLoader.loadPurpleContext;
import static com.hartwig.hmftools.neo.score.DataLoader.loadRnaNeoData;
import static com.hartwig.hmftools.neo.score.DataLoader.loadSomaticVariants;
import static com.hartwig.hmftools.neo.score.NeoScorerConfig.RNA_SAMPLE_ID_SUFFIX;
import static com.hartwig.hmftools.neo.score.TpmCalculator.DEFAULT_PEPTIDE_LENGTH_RANGE;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.isofox.TranscriptExpressionLoader;
import com.hartwig.hmftools.common.neo.RnaNeoEpitope;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.rna.TranscriptExpressionFile;
import com.hartwig.hmftools.neo.bind.BindData;

public class NeoScorerTask implements Callable
{
    private final int mThreadId;
    private final List<SampleData> mSamples;

    private final NeoScorerConfig mConfig;
    private final ReferenceData mReferenceData;
    private final NeoDataWriter mWriters;

    public NeoScorerTask(
            final int threadId, final NeoScorerConfig config, final ReferenceData referenceData, final NeoDataWriter writers)
    {
        mThreadId = threadId;
        mSamples = Lists.newArrayList();
        mConfig = config;
        mReferenceData = referenceData;
        mWriters = writers;
    }

    public void addSample(final SampleData sampleData) { mSamples.add(sampleData); }

    @Override
    public Long call()
    {
        if(mSamples.size() > 1)
        {
            NE_LOGGER.info("{}: processing {} samples", mThreadId, mSamples.size());
        }

        int sampleIndex = 0;

        try
        {
            for(; sampleIndex < mSamples.size(); ++sampleIndex)
            {
                processSample(mSamples.get(sampleIndex));

                if(sampleIndex > 0 && (sampleIndex % 10) == 0)
                {
                    NE_LOGGER.info("{}: processed {} samples of {}", mThreadId, sampleIndex, mSamples.size());
                }
            }
        }
        catch(Exception e)
        {
            NE_LOGGER.error("sample({}) failed processing", mSamples.get(sampleIndex).TumorId, e.toString());
            e.printStackTrace();
            System.exit(1);
        }

        if(mSamples.size() > 1)
        {
            NE_LOGGER.info("{}: processing complete", mThreadId, mSamples.size());
        }

        return (long)1;
    }

    public void processSample(final SampleData sample)
    {
        String sampleId = sample.TumorId;

        List<NeoEpitopeData> neoDataList = loadNeoEpitopes(sampleId, mConfig.NeoDir);

        List<AlleleCoverage> alleleCoverages = loadAlleleCoverage(sampleId, mConfig.LilacDir);

        if(neoDataList == null || alleleCoverages == null)
            System.exit(1);

        NE_LOGGER.info("sample({}) processing {} neoepitopes", sampleId, neoDataList.size());

        Set<String> uniqueAlleles = alleleCoverages.stream().map(x -> x.Allele).collect(Collectors.toSet());

        PurityContext purityContext = loadPurpleContext(mConfig.PurpleDir, sampleId);
        double samplePloidy = purityContext.bestFit().ploidy();

        TpmSource tpmSource = null;

        Map<String,Double> sampleTPMs = Maps.newHashMap();

        if(sample.HasRna)
        {
            if(mReferenceData.TranscriptExpression == null)
            {
                NE_LOGGER.debug("sample({}) loading transcript expression", sampleId);

                try
                {
                    String tpmFileSampleId = sample.RnaSampleId != null ? sample.RnaSampleId : sample.TumorId;
                    sampleTPMs.putAll(TranscriptExpressionLoader.loadTranscriptExpression(mConfig.IsofoxDir, tpmFileSampleId));
                }
                catch(Exception e)
                {
                    NE_LOGGER.error("failed to load sample({}) transcript expression", sampleId, e.toString());
                    System.exit(1);
                }
            }
            else
            {
                if(!mReferenceData.TranscriptExpression.hasSampleId(sampleId))
                {
                    NE_LOGGER.error("sample({}) missing from transcript expression matrix", sampleId);
                    System.exit(1);
                }
            }

            tpmSource = TpmSource.SAMPLE;
        }
        else
        {
            if(mReferenceData.TpmMedians.hasCancerType(sample.CancerType))
                tpmSource = TpmSource.CANCER_TYPE;
            else
                tpmSource = TpmSource.COHORT;
        }

        List<RnaNeoEpitope> rnaNeoDataList = loadRnaNeoData(sample, mConfig.IsofoxDir);

        List<SomaticVariant> somaticVariants = Lists.newArrayList();

        if(!mConfig.RnaSomaticVcf.isEmpty())
        {
            List<NeoEpitopeData> pointNeos = neoDataList.stream().filter(x -> x.VariantType.isPointMutation()).collect(Collectors.toList());
            somaticVariants = loadSomaticVariants(sample, mConfig.RnaSomaticVcf, pointNeos);
        }

        if(sample.HasRna && (rnaNeoDataList == null || somaticVariants == null))
        {
            NE_LOGGER.error("sample({}) missing required RNA: fusions({}) variants({})",
                    sampleId, rnaNeoDataList == null ? "missing" : "present", somaticVariants == null ? "missing" : "present");
            System.exit(1);
        }

        // set TPM and RNA fragment & depth as available
        for(NeoEpitopeData neoData : neoDataList)
        {
            // set sample and cohort TPM values
            neoData.setExpressionData(sample, sampleTPMs, mReferenceData.TranscriptExpression, mReferenceData.TpmMedians);

            // set RNA fragment counts for the specific variant or neoepitope
            neoData.setFusionRnaSupport(rnaNeoDataList);
            neoData.setMutationRnaSupport(somaticVariants);
        }

        TpmCalculator tpmCalculator = new TpmCalculator(FLANK_AA_COUNT, DEFAULT_PEPTIDE_LENGTH_RANGE);

        tpmCalculator.compute(sampleId, neoDataList, samplePloidy);

        // build out results per allele and score them
        int scoreCount = 0;

        for(NeoEpitopeData neoData : neoDataList)
        {
            int i = 0;
            while(i < neoData.peptides().size())
            {
                PeptideScoreData peptideScoreData = neoData.peptides().get(i);

                // check for a wild-type match
                if(neoData.Transcripts[FS_UP].stream().anyMatch(x -> mReferenceData.peptideMatchesWildtype(peptideScoreData.Peptide, x)))
                {
                    neoData.peptides().remove(i);
                    continue;
                }

                if(peptideScoreData.alleleScoreData().isEmpty())
                {
                    uniqueAlleles.forEach(x -> peptideScoreData.addAllele(x));

                    for(BindData bindData : peptideScoreData.alleleScoreData())
                    {
                        mReferenceData.PeptideScorer.calcScoreData(bindData);
                        ++scoreCount;
                    }
                }

                ++i;
            }
        }

        NE_LOGGER.debug("sample({}) neoepitopes({}) scored {} allele-peptides",
                sampleId, neoDataList.size(), scoreCount);

        if(mConfig.WriteTypes.contains(OutputType.ALLELE_PEPTIDE))
        {
            for(NeoEpitopeData neoData : neoDataList)
            {
                mWriters.writePeptideData(sampleId, neoData, alleleCoverages);
            }
        }

        if(mConfig.WriteTypes.contains(OutputType.NEOEPITOPE))
        {
            TpmSource sampleTmpSource = tpmSource;
            neoDataList.forEach(x -> mWriters.writeNeoData(sampleId, sampleTmpSource, x));
        }
    }
}
