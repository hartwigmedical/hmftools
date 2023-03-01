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
    private final String mSampleId;
    private final NeoScorerConfig mConfig;
    private final BindScorer mScorer;
    private final RnaExpressionMatrix mTransExpression;
    private final NeoDataWriter mWriters;

    private static final int[] PEPTIDE_LENGTH_RANGE = new int[] { MIN_PEPTIDE_LENGTH, REF_PEPTIDE_LENGTH };

    public NeoScorerTask(
            final String sampleId, final NeoScorerConfig config, final BindScorer scorer,
            @Nullable final RnaExpressionMatrix transExpression, final NeoDataWriter writers)
    {
        mSampleId = sampleId;
        mConfig = config;
        mScorer = scorer;
        mTransExpression = transExpression;
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
            NE_LOGGER.error("sample({}) failed processing", mSampleId, e.toString());
            e.printStackTrace();
            System.exit(1);
        }

        return (long)1;
    }

    public void processSample() throws Exception
    {
        Map<Integer,NeoEpitopeData> neoEpitopeMap = loadNeoEpitopes(mSampleId, mConfig.NeoDataDir);

        List<AlleleCoverage> alleleCoverages = loadAlleleCoverage(mSampleId, mConfig.LilacDataDir);

        List<RnaNeoEpitope> rnaNeoDataList = loadRnaNeoData(mSampleId, mConfig.IsofoxDataDir);

        Map<String,Double> sampleTPMs = Maps.newHashMap();

        if(mTransExpression == null || !mTransExpression.hasSampleId(mSampleId))
        {
            sampleTPMs.putAll(TranscriptExpressionLoader.loadTranscriptExpression(mConfig.IsofoxDataDir, mSampleId));
        }

        List<SomaticVariant> somaticVariants = Lists.newArrayList();

        if(!mConfig.RnaSomaticVcf.isEmpty())
        {
            String rnaSampleId = mSampleId + RNA_SAMPLE_APPEND_SUFFIX;
            List<NeoEpitopeData> pointNeos = neoEpitopeMap.values().stream().filter(x -> x.VariantType.isPointMutation()).collect(Collectors.toList());
            somaticVariants.addAll(loadSomaticVariants(mSampleId, rnaSampleId, mConfig.RnaSomaticVcf, pointNeos));
        }

        // set TPM and RNA fragment & depth as available
        for(NeoEpitopeData neoData : neoEpitopeMap.values())
        {
            neoData.setExpressionData(mSampleId, sampleTPMs, mTransExpression);
            neoData.setFusionRnaSupport(rnaNeoDataList);
            neoData.setMutationRnaSupport(somaticVariants);
        }

        double tpmNormalisationFactor = calculateTpmNormalisation(neoEpitopeMap.values());

        Map<Integer,NeoPredictionData> neoPredictionsMap = Maps.newHashMap();

        int peptideAlleleCount = 0;

        for(NeoEpitopeData neoData : neoEpitopeMap.values())
        {
            // derive the set of peptides per allele from the novel amino acids
            NeoPredictionData predData = produceAllelePeptides(neoData, alleleCoverages);

            neoPredictionsMap.put(neoData.Id, predData);

            for(List<BindData> bindDataList : predData.getPeptidePredictions().values())
            {
                for(BindData bindData : bindDataList)
                {
                    mScorer.calcScoreData(bindData);
                }
            }

            peptideAlleleCount += predData.getPeptidePredictions().values().stream().mapToInt(x -> x.size()).sum();
        }

        NE_LOGGER.debug("sample({}) neoepitopes({}) scored {} allele-peptides",
                mSampleId, neoEpitopeMap.size(), peptideAlleleCount);

        if(mConfig.WriteTypes.contains(OutputType.ALLELE_PEPTIDE))
        {
            for(AlleleCoverage alleleCoverage : alleleCoverages)
            {
                for(NeoEpitopeData neoData : neoEpitopeMap.values())
                {
                    NeoPredictionData predData = neoPredictionsMap.get(neoData.Id);

                    mWriters.writePeptideData(mSampleId, neoData, predData, alleleCoverage);
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
