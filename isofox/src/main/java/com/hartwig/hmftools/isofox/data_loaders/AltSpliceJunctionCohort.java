package com.hartwig.hmftools.isofox.data_loaders;

import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.data_loaders.DataLoadType.ALT_SPLICE_JUNCTIONS;
import static com.hartwig.hmftools.isofox.data_loaders.DataLoaderConfig.formSampleFilenames;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.stats.FisherExactTest;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunction;

public class AltSpliceJunctionCohort
{
    private final DataLoaderConfig mConfig;

    // map of chromosomes to a map of genes to a list of alternate splice junctions
    private final Map<String, Map<String,List<AltSpliceJunction>>> mAltSpliceJunctions;

    private final FisherExactTest mFisherET;

    public AltSpliceJunctionCohort(final DataLoaderConfig config)
    {
        mConfig = config;
        mAltSpliceJunctions = Maps.newHashMap();
        mFisherET = new FisherExactTest();
    }

    public void processAltSpliceJunctions()
    {
        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, ALT_SPLICE_JUNCTIONS, filenames))
            return;

        // load each sample's alt SJs and consolidate into a single list
        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path altSJFile = filenames.get(i);

            final List<AltSpliceJunction> altSJs = AltSpliceJunction.loadFile(altSJFile);

            ISF_LOGGER.debug("{}: sample({}) loaded {} alt-SJ records", i, sampleId, altSJs.size());

            altSJs.forEach(x -> addAltSpliceJunction(x, sampleId));
        }

        // write a report for any re-occurring alt SJ
        writeReoccurringAltSpliceJunctions();
    }

    private void addAltSpliceJunction(final AltSpliceJunction altSJ, final String sampleId)
    {
        if(!mConfig.RestrictedGeneIds.isEmpty() && !mConfig.RestrictedGeneIds.contains(altSJ.getGeneId()))
            return;

        if(!mConfig.ExcludedGeneIds.isEmpty() && mConfig.ExcludedGeneIds.contains(altSJ.getGeneId()))
            return;

        Map<String,List<AltSpliceJunction>> chrSJs = mAltSpliceJunctions.get(altSJ.Chromosome);

        if(chrSJs == null)
        {
            chrSJs = Maps.newHashMap();
            mAltSpliceJunctions.put(altSJ.Chromosome, chrSJs);
        }

        List<AltSpliceJunction> geneList = chrSJs.get(altSJ.getGeneId());

        if(geneList == null)
        {
            geneList = Lists.newArrayList();
            chrSJs.put(altSJ.getGeneId(), geneList);
        }

        AltSpliceJunction existingAltSJ = geneList.stream().filter(x -> x.matches(altSJ)).findFirst().orElse(null);

        if(existingAltSJ != null)
        {
            existingAltSJ.getSampleIds().add(sampleId);
            existingAltSJ.addFragmentCount(altSJ.getFragmentCount());
            existingAltSJ.addPositionCount(SE_START, altSJ.getPositionCount(SE_START));
            existingAltSJ.addPositionCount(SE_END, altSJ.getPositionCount(SE_END));
            return;
        }

        altSJ.getSampleIds().add(sampleId);
        geneList.add(altSJ);
    }

    private void writeReoccurringAltSpliceJunctions()
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("alt_sj_reoccurring.csv");
            final BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("GeneId,Chromosome,Type,SjStart,SjEnd");
            writer.write(",StartContext,EndContext,BaseMotif,AvgFragCount,MaxFragCount,AvgStartDepth,AvgEndDepth");
            writer.write(",SampleCount,FetProb,ExpVal,AltSJAndMutation,AltSJNoMutation");
            writer.write(",WithMutation,NoAltSJWithMutation,NoAltSJNoMutation");
            writer.write(",SampleIds");
            writer.newLine();

            int totalSampleCount = mConfig.SampleData.SampleIds.size();
            mFisherET.initialise(totalSampleCount);

            for(Map.Entry<String,Map<String,List<AltSpliceJunction>>> chrEntry : mAltSpliceJunctions.entrySet())
            {
                final String chromosome = chrEntry.getKey();
                final Map<String,List<AltSpliceJunction>> geneMap = chrEntry.getValue();

                for(Map.Entry<String,List<AltSpliceJunction>> geneEntry : geneMap.entrySet())
                {
                    final String geneId = geneEntry.getKey();

                    for (AltSpliceJunction altSJ : geneEntry.getValue())
                    {
                        int scWithAltSJ = altSJ.getSampleIds().size();

                        if(scWithAltSJ < mConfig.AltSJMinSampleThreshold)
                            continue;

                        writer.write(String.format("%s,%s,%s,%d,%d",
                                geneId, chromosome, altSJ.type(),
                                altSJ.SpliceJunction[SE_START], altSJ.SpliceJunction[SE_END]));

                        writer.write(String.format(",%s,%s,%s,%.0f,%d,%.0f,%.0f",
                                altSJ.RegionContexts[SE_START], altSJ.RegionContexts[SE_END], altSJ.getDonorAcceptorBases(),
                                altSJ.getFragmentCount()/(double)scWithAltSJ, altSJ.getMaxFragmentCount(),
                                altSJ.getPositionCount(SE_START)/(double)scWithAltSJ,
                                altSJ.getPositionCount(SE_END)/(double)scWithAltSJ));

                        int scWithAltSJWithMutation = mConfig.SampleData.sampleCountWithMutation(altSJ.getSampleIds());
                        int scWithAltSJNoMutation = scWithAltSJ - scWithAltSJWithMutation;
                        int scWithMutation = mConfig.SampleData.SampleMutationType.size();
                        int scNoAltSJWithMutation = scWithMutation - scWithAltSJWithMutation;
                        int scNoAltSJNoMutation = totalSampleCount - scWithAltSJWithMutation - scWithAltSJNoMutation - scNoAltSJWithMutation;

                        double expectedVal  = scWithMutation * scWithAltSJ / (double)totalSampleCount;
                        double fisherProb = mFisherET.calc(scWithAltSJWithMutation, scNoAltSJWithMutation, scWithAltSJNoMutation, scNoAltSJNoMutation, expectedVal);

                        writer.write(String.format(",%d,%4.3e,%.1f,%d,%d,%d,%d,%d",
                                scWithAltSJ, fisherProb, expectedVal,
                                scWithAltSJWithMutation, scWithAltSJNoMutation, scWithMutation, scNoAltSJWithMutation, scNoAltSJNoMutation));

                        writer.write(String.format(",%s", appendStrList(altSJ.getSampleIds(), ';')));

                        writer.newLine();
                    }
                }
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write cohort summary data file: {}", e.toString());
        }

    }
}
