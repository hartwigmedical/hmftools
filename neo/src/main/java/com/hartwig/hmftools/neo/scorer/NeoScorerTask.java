package com.hartwig.hmftools.neo.scorer;

import static com.hartwig.hmftools.common.codon.AminoAcidRna.AA_SELENOCYSTEINE;
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

    private static final int[] PEPTIDE_LENGTH_RANGE = new int[] { MIN_PEPTIDE_LENGTH, REF_PEPTIDE_LENGTH };

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
            NE_LOGGER.error("sample({}) failed processing", mSampleData, e.toString());
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

        List<RnaNeoEpitope> rnaNeoDataList = loadRnaNeoData(sampleId, mConfig.IsofoxDataDir);

        Map<String,Double> sampleTPMs = Maps.newHashMap();

        if(mTransExpression == null || !mTransExpression.hasSampleId(sampleId))
        {
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

        double tpmNormalisationFactor = calculateTpmNormalisation(neoDataList);

        // derive the set of peptides per allele from the novel amino acids
        Map<Integer,NeoPredictionData> neoPredictionsMap = Maps.newHashMap();

        int peptideAlleleCount = 0;

        Map<String,List<NeoEpitopeData>> variantNeoMap = Maps.newHashMap();

        for(NeoEpitopeData neoData : neoDataList)
        {
            List<NeoEpitopeData> variantNeos = variantNeoMap.get(neoData.VariantInfo);

            if(variantNeos == null)
                variantNeoMap.put(neoData.VariantInfo, Lists.newArrayList(neoData));
            else
                variantNeos.add(neoData);
        }

        for(List<NeoEpitopeData> variantNeos : variantNeoMap.values())
        {

            for(NeoEpitopeData neoData : neoDataList)
            {
                NeoPredictionData predData = produceAllelePeptides(neoData, alleleCoverages);

                neoPredictionsMap.put(neoData.Id, predData);
                peptideAlleleCount += predData.getPeptidePredictions().values().stream().mapToInt(x -> x.size()).sum();

                // sum TPM across all transcripts

            }
        }

        NE_LOGGER.debug("sample({}) neoepitopes({}) derived {} allele-peptides",
                sampleId, neoDataList.size(), peptideAlleleCount);

        // calculate TPM per allele and peptide



        // score the peptides
        for(NeoPredictionData predData : neoPredictionsMap.values())
        {
            for(List<BindData> bindDataList : predData.getPeptidePredictions().values())
            {
                for(BindData bindData : bindDataList)
                {
                    mScorer.calcScoreData(bindData);
                }
            }
        }

        NE_LOGGER.debug("sample({}) neoepitopes({}) scored {} allele-peptides",
                sampleId, neoDataList.size(), peptideAlleleCount);

        if(mConfig.WriteTypes.contains(OutputType.ALLELE_PEPTIDE))
        {
            for(AlleleCoverage alleleCoverage : alleleCoverages)
            {
                for(NeoEpitopeData neoData : neoDataList)
                {
                    NeoPredictionData predData = neoPredictionsMap.get(neoData.Id);

                    mWriters.writePeptideData(sampleId, neoData, predData, alleleCoverage);
                }
            }
        }
    }

    private double calculateTpmNormalisation(final Collection<NeoEpitopeData> neoDataList)
    {
        double tpmUpTotal = 0;
        double baseDepthTotal = 0;
        int dataCount = 0;

        for(NeoEpitopeData neoData : neoDataList)
        {
            if(neoData.VariantType != NeoEpitopeType.MISSENSE)
                continue;

            if(neoData.RnaData.hasCoverage() && neoData.RnaData.hasExpression())
            {
                tpmUpTotal += neoData.RnaData.transExpression()[FS_UP];
                baseDepthTotal += neoData.RnaData.averageBaseDepth();
                ++dataCount;
            }
        }

        double tpmNormalisationFactor = dataCount > 0 ? tpmUpTotal / baseDepthTotal : 0;
        return tpmNormalisationFactor;
    }

    private NeoPredictionData produceAllelePeptides(final NeoEpitopeData neoData, final List<AlleleCoverage> alleleCoverages)
    {
        NeoPredictionData neoPredData = new NeoPredictionData(neoData.Id);

        final List<PeptideData> peptides = EpitopeUtils.generatePeptides(
                neoData.UpAminoAcids, neoData.NovelAminoAcids, neoData.DownAminoAcids, PEPTIDE_LENGTH_RANGE, FLANK_AA_COUNT);

        Set<String> uniqueAlleles = Sets.newHashSet();

        for(AlleleCoverage allele : alleleCoverages)
        {
            if(uniqueAlleles.contains(allele.Allele))
                continue;

            uniqueAlleles.add(allele.Allele);

            List<BindData> bindDataList = Lists.newArrayList();
            neoPredData.getPeptidePredictions().put(allele.Allele, bindDataList);

            for(PeptideData peptideData : peptides)
            {
                if(peptideData.Peptide.contains(AA_SELENOCYSTEINE))
                    continue;

                BindData bindData = new BindData(
                        allele.Allele, peptideData.Peptide, "", peptideData.UpFlank, peptideData.DownFlank);

                bindData.setTPM(neoData.RnaData.getTPM());

                bindDataList.add(bindData);
            }
        }

        return neoPredData;
    }
}
