package com.hartwig.hmftools.neo.scorer;

import static com.hartwig.hmftools.common.codon.AminoAcidRna.AA_SELENOCYSTEINE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.MIN_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.FlankCounts.FLANK_AA_COUNT;
import static com.hartwig.hmftools.neo.scorer.DataLoader.loadAlleleCoverage;
import static com.hartwig.hmftools.neo.scorer.DataLoader.loadNeoEpitopes;
import static com.hartwig.hmftools.neo.scorer.DataLoader.loadRnaNeoData;
import static com.hartwig.hmftools.neo.scorer.DataLoader.loadSomaticVariants;
import static com.hartwig.hmftools.neo.scorer.NeoScorerConfig.RNA_SAMPLE_APPEND_SUFFIX;
import static com.hartwig.hmftools.neo.scorer.TpmCalculator.DEFAULT_PEPTIDE_LENGTH_RANGE;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.isofox.TranscriptExpressionLoader;
import com.hartwig.hmftools.common.neo.NeoEpitopeType;
import com.hartwig.hmftools.common.neo.RnaNeoEpitope;
import com.hartwig.hmftools.common.rna.RnaExpressionMatrix;
import com.hartwig.hmftools.neo.bind.BindData;
import com.hartwig.hmftools.neo.bind.BindScorer;
import com.hartwig.hmftools.neo.epitope.EpitopeUtils;
import com.hartwig.hmftools.neo.PeptideData;

import org.jetbrains.annotations.Nullable;

public class NeoScorerTask implements Callable
{
    private final SampleData mSampleData;
    private final NeoScorerConfig mConfig;
    private final BindScorer mScorer;
    private final RnaExpressionMatrix mTransExpression;
    private final TpmMediansCache mTpmMediansCache;
    private final NeoDataWriter mWriters;

    public NeoScorerTask(
            final SampleData sampleData, final NeoScorerConfig config, final BindScorer scorer,
            @Nullable final RnaExpressionMatrix transExpression, final TpmMediansCache tpmMediansCache, final NeoDataWriter writers)
    {
        mSampleData = sampleData;
        mConfig = config;
        mScorer = scorer;
        mTransExpression = transExpression;
        mTpmMediansCache = tpmMediansCache;
        mWriters = writers;
    }

    @Override
    public Long call()
    {
        try
        {
            processSample();
        }
        catch(Exception e)
        {
            NE_LOGGER.error("sample({}) failed processing", mSampleData.Id, e.toString());
            e.printStackTrace();
            System.exit(1);
        }

        return (long)1;
    }

    public void processSample() throws Exception
    {
        String sampleId = mSampleData.Id;

        List<NeoEpitopeData> neoDataList = loadNeoEpitopes(sampleId, mConfig.NeoDataDir);

        List<AlleleCoverage> alleleCoverages = loadAlleleCoverage(sampleId, mConfig.LilacDataDir);
        Set<String> uniqueAlleles = alleleCoverages.stream().map(x -> x.Allele).collect(Collectors.toSet());

        List<RnaNeoEpitope> rnaNeoDataList = loadRnaNeoData(sampleId, mConfig.IsofoxDataDir);

        Map<String,Double> sampleTPMs = Maps.newHashMap();

        if(mTransExpression == null || !mTransExpression.hasSampleId(sampleId))
        {
            NE_LOGGER.debug("sample({}) loading transcript expression", sampleId);
            sampleTPMs.putAll(TranscriptExpressionLoader.loadTranscriptExpression(mConfig.IsofoxDataDir, sampleId));
        }

        List<SomaticVariant> somaticVariants = Lists.newArrayList();

        if(!mConfig.RnaSomaticVcf.isEmpty())
        {
            String rnaSampleId = sampleId + RNA_SAMPLE_APPEND_SUFFIX;
            List<NeoEpitopeData> pointNeos = neoDataList.stream().filter(x -> x.VariantType.isPointMutation()).collect(Collectors.toList());
            somaticVariants.addAll(loadSomaticVariants(sampleId, rnaSampleId, mConfig.RnaSomaticVcf, pointNeos));
        }

        // set TPM and RNA fragment & depth as available
        for(NeoEpitopeData neoData : neoDataList)
        {
            // set sample and cohort TPM values
            neoData.setExpressionData(mSampleData, sampleTPMs, mTransExpression, mTpmMediansCache);

            // set RNA fragment counts for the specific variant or neoepitope
            neoData.setFusionRnaSupport(rnaNeoDataList);
            neoData.setMutationRnaSupport(somaticVariants);
        }

        TpmCalculator tpmCalculator = new TpmCalculator(FLANK_AA_COUNT, DEFAULT_PEPTIDE_LENGTH_RANGE);

        tpmCalculator.compute(sampleId, neoDataList);

        // build out results per allele and score them
        int scoreCount = 0;

        for(NeoEpitopeData neoData : neoDataList)
        {
            for(PeptideScoreData peptideScoreData : neoData.peptides())
            {
                if(!peptideScoreData.alleleScoreData().isEmpty())
                    continue;

                uniqueAlleles.forEach(x -> peptideScoreData.addAllele(x));

                for(BindData bindData : peptideScoreData.alleleScoreData())
                {
                    mScorer.calcScoreData(bindData);
                    ++scoreCount;
                }
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
            neoDataList.forEach(x -> mWriters.writeNeoData(sampleId, x));
        }
    }
}
