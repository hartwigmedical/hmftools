package com.hartwig.hmftools.neo.missense;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.neo.NeoCommon.APP_NAME;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindScorer.INVALID_CALC;
import static com.hartwig.hmftools.neo.bind.ScoreConfig.SCORE_FILE_DIR;
import static com.hartwig.hmftools.neo.missense.MissenseConfig.registerConfig;
import static com.hartwig.hmftools.neo.score.NeoRnaData.NO_TPM_VALUE;
import static com.hartwig.hmftools.neo.score.NeoScorerConfig.COHORT_TPM_MEDIANS_FILE;
import static com.hartwig.hmftools.neo.score.TpmMediansCache.PAN_CANCER_VALUE;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.neo.bind.BindData;
import com.hartwig.hmftools.neo.bind.BindScorer;
import com.hartwig.hmftools.neo.bind.ScoreConfig;
import com.hartwig.hmftools.neo.score.TpmMediansCache;

import org.jetbrains.annotations.NotNull;

public class MissensePeptideScorer
{
    private final EnsemblDataCache mEnsemblDataCache;
    private final RefGenomeInterface mRefGenome;
    private final MissenseConfig mConfig;

    private final BindScorer mPeptideScorer;
    private final MissenseWriter mWriter;
    private final TpmMediansCache mCohortTpmMedians;

    private final List<String> mAlleleList;

    public MissensePeptideScorer(final ConfigBuilder configBuilder)
    {
        mEnsemblDataCache = new EnsemblDataCache(configBuilder);
        mEnsemblDataCache.setRequiredData(true, false, false, true);

        mRefGenome = loadRefGenome(configBuilder.getValue(REF_GENOME));

        mConfig = new MissenseConfig(configBuilder);

        mPeptideScorer = configBuilder.hasValue(SCORE_FILE_DIR) ? new BindScorer(new ScoreConfig(configBuilder)) : null;

        mAlleleList = Lists.newArrayList();

        mCohortTpmMedians = new TpmMediansCache(configBuilder.getValue(COHORT_TPM_MEDIANS_FILE));

        mWriter = new MissenseWriter(mConfig, mPeptideScorer != null);
    }

    public void run()
    {
        if(mConfig.GeneIds.isEmpty())
        {
            NE_LOGGER.error("no gene IDs specified");
            System.exit(1);
        }

        if(mPeptideScorer != null && !mPeptideScorer.loadScoringData())
        {
            NE_LOGGER.error("bind scoring data load failed");
            System.exit(1);
        }

        long startTimeMs = System.currentTimeMillis();

        mEnsemblDataCache.load(true);
        mEnsemblDataCache.loadTranscriptData(mConfig.GeneIds);

        if(!mConfig.Alleles.isEmpty())
            mAlleleList.addAll(mConfig.Alleles);
        else if(mPeptideScorer != null)
            mPeptideScorer.getScoringAlleles().forEach(x -> mAlleleList.add(x));

        NE_LOGGER.info("generating missense peptides for {} genes", mConfig.GeneIds.size());

        if(mPeptideScorer != null)
        {
            NE_LOGGER.info("scoring peptides for {} alleles", mAlleleList.size());
        }

        List<GeneTask> geneTasks = Lists.newArrayList();

        if(mConfig.Threads > 1)
        {
            boolean allocateByGene = mConfig.GeneIds.size() >= mAlleleList.size();
            int maxItems = allocateByGene ? mConfig.GeneIds.size() : mAlleleList.size();

            for(int i = 0; i < min(maxItems, mConfig.Threads); ++i)
            {
                geneTasks.add(new GeneTask(i));
            }

            int taskIndex = 0;

            if(allocateByGene)
            {
                geneTasks.forEach(x -> x.Alleles.addAll(mAlleleList));

                for(String geneId : mConfig.GeneIds)
                {
                    if(taskIndex >= geneTasks.size())
                        taskIndex = 0;

                    geneTasks.get(taskIndex).GeneIds.add(geneId);
                    ++taskIndex;
                }
            }
            else
            {
                geneTasks.forEach(x -> x.GeneIds.addAll(mConfig.GeneIds));

                for(String allele : mAlleleList)
                {
                    if(taskIndex >= geneTasks.size())
                        taskIndex = 0;

                    geneTasks.get(taskIndex).Alleles.add(allele);
                    ++taskIndex;
                }
            }

            final List<Callable<Void>> callableList = geneTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, mConfig.Threads);
        }
        else
        {
            GeneTask geneTask = new GeneTask(0);
            geneTask.GeneIds.addAll(mConfig.GeneIds);
            geneTask.Alleles.addAll(mAlleleList);
            geneTask.call();
        }

        mWriter.close();

        NE_LOGGER.info("missense peptide generation complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private class GeneTask implements Callable<Void>
    {
        public final List<String> GeneIds;
        public final List<String> Alleles;
        private final int mTaskId;
        private final MissenseCalcs mMissenseCalcs;

        public GeneTask(int taskId)
        {
            GeneIds = Lists.newArrayList();
            Alleles = Lists.newArrayList();
            mTaskId = taskId;
            mMissenseCalcs = new MissenseCalcs(mConfig, mRefGenome);
        }

        @Override
        public Void call()
        {
            int geneCount = 0;

            NE_LOGGER.info("{}: processing {} genes for {} alleles", mTaskId, GeneIds.size(), Alleles.size());

            for(String geneId : GeneIds)
            {
                GeneData geneData = mEnsemblDataCache.getGeneDataById(geneId);

                List<TranscriptData> transDataList = mEnsemblDataCache.getTranscripts(geneData.GeneId);

                if(transDataList == null)
                    continue;

                for(TranscriptData transData : transDataList)
                {
                    if(!transData.IsCanonical)
                        continue;

                    mMissenseCalcs.clear();

                    mMissenseCalcs.processTranscript(geneData, transData);

                    NE_LOGGER.debug("gene({}:{}) scoring {} missense peptides",
                            geneData.GeneId, geneData.GeneName, mMissenseCalcs.peptideData().size());

                    if(mPeptideScorer != null)
                    {
                        scorePeptides(mMissenseCalcs.peptideData(), Alleles);
                    }
                    else
                    {
                        mWriter.writePeptideData(mMissenseCalcs.peptideData(), null);
                    }
                }

                ++geneCount;

                if(geneCount > 0 && (geneCount % 10) == 0)
                {
                    NE_LOGGER.info("{}: processed {} genes", mTaskId, geneCount);
                }
            }

            return null;
        }
    }

    private void scorePeptides(final List<MissensePeptide> peptideData, final List<String> alleles)
    {
        if(peptideData.isEmpty())
            return;

        // NOTE: peptides are processed gene by gene, so all have the same transcript
        double[] tpmValues = mCohortTpmMedians.getTranscriptTpm(peptideData.get(0).TransName, null);

        for(String allele : alleles)
        {
            NE_LOGGER.debug("gene({}) allele({}) calculating binding scores for {} peptides",
                    peptideData.get(0).GeneName, allele, peptideData.size());

            final List<MissensePeptide> peptideDataList = Lists.newArrayList();
            final List<BindData> bindDataList = Lists.newArrayList();

            for(MissensePeptide misensePeptide : peptideData)
            {
                BindData bindData = new BindData(allele, misensePeptide.Peptide, "", misensePeptide.UpFlank, misensePeptide.DownFlank);

                if(tpmValues[PAN_CANCER_VALUE] != NO_TPM_VALUE)
                    bindData.setTPM(tpmValues[PAN_CANCER_VALUE]);

                mPeptideScorer.calcScoreData(bindData);

                if(passesRankThreshold(bindData.likelihoodRank()))
                {
                    peptideDataList.add(misensePeptide);
                    bindDataList.add(bindData);
                }
            }

            mWriter.writePeptideData(peptideDataList, bindDataList);
        }
    }

    private boolean passesRankThreshold(double value)
    {
        return value != INVALID_CALC && (mConfig.LikelihoodCutoff == 0 || mConfig.LikelihoodCutoff > 0 && value <= mConfig.LikelihoodCutoff);
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        MissensePeptideScorer missensePeptideScorer = new MissensePeptideScorer(configBuilder);
        missensePeptideScorer.run();
    }
}
