package com.hartwig.hmftools.neo.cohort;

import static com.hartwig.hmftools.common.codon.AminoAcidRna.AA_SELENOCYSTEINE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.rna.RnaExpressionMatrix.INVALID_EXP;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.MIN_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.FlankCounts.FLANK_BASE_COUNT;
import static com.hartwig.hmftools.neo.cohort.DataLoader.loadAlleleCoverage;
import static com.hartwig.hmftools.neo.cohort.DataLoader.loadNeoEpitopes;
import static com.hartwig.hmftools.neo.cohort.DataLoader.loadRnaNeoData;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.neo.RnaNeoEpitope;
import com.hartwig.hmftools.common.rna.RnaExpressionMatrix;
import com.hartwig.hmftools.neo.bind.BindData;
import com.hartwig.hmftools.neo.bind.BindScorer;
import com.hartwig.hmftools.neo.epitope.EpitopeUtils;
import com.hartwig.hmftools.neo.PeptideData;

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
            final RnaExpressionMatrix transExpression, final NeoDataWriter writers)
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
        processSample();
        return (long)1;
    }

    public void processSample()
    {
        Map<Integer,NeoEpitopeData> neoEpitopeMap = loadNeoEpitopes(mSampleId, mConfig.NeoDataDir);

        List<AlleleCoverage> alleleCoverages = loadAlleleCoverage(mSampleId, mConfig.LilacDataDir);

        List<RnaNeoEpitope> rnaNeoDataList = loadRnaNeoData(mSampleId, mConfig.IsofoxDataDir);

        Map<Integer,NeoPredictionData> neoPredictionsMap = Maps.newHashMap();

        int peptideAlleleCount = 0;

        for(Map.Entry<Integer,NeoEpitopeData> entry : neoEpitopeMap.entrySet())
        {
            NeoEpitopeData neoData = entry.getValue();

            if(mTransExpression.hasSampleId(mSampleId))
            {
                for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
                {
                    neoData.TransExpression[fs] = 0;

                    for(String transName : neoData.Transcripts[fs])
                    {
                        double expression = mTransExpression.getExpression(transName, mSampleId);

                        // distinguish non-existent expression vs zero TPM
                        if(expression != INVALID_EXP)
                            neoData.TransExpression[fs] += expression;
                    }
                }
            }

            RnaNeoEpitope rnaNeoData = rnaNeoDataList.stream()
                    .filter(x -> x.Id == neoData.Id && x.VariantInfo.equals(neoData.VariantInfo)).findFirst().orElse(null);

            if(rnaNeoData != null)
            {
                neoData.RnaNovelFragments = rnaNeoData.FragmentCount;
                neoData.RnaBaseDepth[FS_UP] = rnaNeoData.BaseDepth[FS_UP];
                neoData.RnaBaseDepth[FS_DOWN] = rnaNeoData.BaseDepth[FS_DOWN];
            }

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

    private NeoPredictionData produceAllelePeptides(final NeoEpitopeData neoData, final List<AlleleCoverage> alleleCoverages)
    {
        NeoPredictionData neoPredData = new NeoPredictionData(neoData.Id);

        final List<PeptideData> peptides = EpitopeUtils.generatePeptides(
                neoData.UpAminoAcids, neoData.NovelAminoAcids, neoData.DownAminoAcids, PEPTIDE_LENGTH_RANGE, FLANK_BASE_COUNT);

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

                bindData.setTPM(neoData.getTPM());

                bindDataList.add(bindData);
            }
        }

        return neoPredData;
    }
}
