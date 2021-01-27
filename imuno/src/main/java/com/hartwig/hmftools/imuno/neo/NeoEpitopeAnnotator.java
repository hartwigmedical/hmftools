package com.hartwig.hmftools.imuno.neo;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_NONSENSE_MED_DECAY;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFusion.DELIMITER;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFusion.NE_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.imuno.common.ImunoCommon.DOWNSTREAM_PRE_GENE_DISTANCE;
import static com.hartwig.hmftools.imuno.common.ImunoCommon.IM_LOGGER;
import static com.hartwig.hmftools.imuno.common.ImunoCommon.LOG_DEBUG;
import static com.hartwig.hmftools.common.neo.AminoAcidConverter.STOP_SYMBOL;
import static com.hartwig.hmftools.imuno.neo.CohortTpmData.CANCER_VALUE;
import static com.hartwig.hmftools.imuno.neo.CohortTpmData.COHORT_VALUE;
import static com.hartwig.hmftools.imuno.neo.NeoConfig.CANCER_TPM_FILE;
import static com.hartwig.hmftools.imuno.neo.NeoConfig.SV_FUSION_FILE;
import static com.hartwig.hmftools.imuno.neo.NeoConfig.GENE_TRANSCRIPTS_DIR;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.neo.NeoEpitopeFile;
import com.hartwig.hmftools.common.neo.NeoEpitopeFusion;
import com.hartwig.hmftools.common.neo.NeoEpitopeType;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;
import org.jooq.Record;
import org.jooq.Result;

public class NeoEpitopeAnnotator
{
    private final NeoConfig mConfig;
    private final Map<String,List<NeoEpitopeFusion>> mSampleFusionMap;

    private final EnsemblDataCache mGeneTransCache;
    private final DatabaseAccess mDbAccess;
    private final CohortTpmData mCohortTpmData;

    private String mCurrentSampleId;
    private BufferedWriter mWriter;

    public NeoEpitopeAnnotator(final CommandLine cmd)
    {
        mConfig = new NeoConfig(cmd);
        mWriter = null;

        mSampleFusionMap = Maps.newHashMap();

        mGeneTransCache = new EnsemblDataCache(cmd.getOptionValue(GENE_TRANSCRIPTS_DIR), RefGenomeVersion.RG_37);
        mGeneTransCache.setRequiredData(true, false, false, false);
        mGeneTransCache.setRestrictedGeneIdList(mConfig.RestrictedGeneIds);
        mGeneTransCache.load(false);
        mGeneTransCache.createGeneNameIdMap();

        mDbAccess = DatabaseAccess.createDatabaseAccess(cmd);
        mCurrentSampleId = "";

        mCohortTpmData = new CohortTpmData(cmd.getOptionValue(CANCER_TPM_FILE));
        loadSvNeoEpitopes(cmd.getOptionValue(SV_FUSION_FILE));
    }

    public void run()
    {
        // check required inputs and config
        for(final String sampleId : mConfig.SampleIds)
        {
            mCurrentSampleId = sampleId;

            final List<NeoEpitopeFusion> fusions = getSvFusions(sampleId);
            final List<PointMutationData> pointMutations = getSomaticVariants(sampleId);

            IM_LOGGER.debug("sample({}) loaded {} fusions and {} point mutations",
                    sampleId, fusions.size(), pointMutations.size());

            addSvFusions(fusions);
            addPointMutations(pointMutations);
        }

        closeBufferedWriter(mWriter);
    }

    private final List<PointMutationData> getSomaticVariants(final String sampleId)
    {
        final List<PointMutationData> pointMutations = Lists.newArrayList();

        final Result<Record> result = mDbAccess.context().select().from(Tables.SOMATICVARIANT)
                .where(Tables.SOMATICVARIANT.SAMPLEID.eq(sampleId))
                .and(Tables.SOMATICVARIANT.FILTER.eq(PASS_FILTER))
                .and(Tables.SOMATICVARIANT.GENE.notEqual(""))
                .fetch();

        for (Record record : result)
        {
            final String gene = record.getValue(Tables.SOMATICVARIANT.GENE);

            if(gene.isEmpty() || mGeneTransCache.getGeneDataByName(gene) == null)
                continue;

            CodingEffect codingEffect = CodingEffect.valueOf(record.getValue(Tables.SOMATICVARIANT.WORSTCODINGEFFECT));

            if(codingEffect != NONSENSE_OR_FRAMESHIFT && codingEffect != MISSENSE)
                continue;

            String chromosome = record.getValue(Tables.SOMATICVARIANT.CHROMOSOME);
            int position = record.getValue(Tables.SOMATICVARIANT.POSITION);
            String ref = record.getValue(Tables.SOMATICVARIANT.REF);
            String alt = record.getValue(Tables.SOMATICVARIANT.ALT);
            double copyNumber = record.getValue(Tables.SOMATICVARIANT.VARIANTCOPYNUMBER);
            Integer localPhaseSet = record.get(Tables.SOMATICVARIANT.LOCALPHASESET);

            pointMutations.add(new PointMutationData(chromosome, position, ref, alt, gene,
                    codingEffect, copyNumber, localPhaseSet != null ? localPhaseSet : -1));
        }

        return pointMutations;
    }

    private final List<NeoEpitopeFusion> getSvFusions(final String sampleId)
    {
        final List<NeoEpitopeFusion> fusions = mSampleFusionMap.get(sampleId);
        mSampleFusionMap.remove(sampleId);
        return fusions != null ? fusions : Lists.newArrayList();
    }

    private void addSvFusions(final List<NeoEpitopeFusion> fusions)
    {
        for(NeoEpitopeFusion fusion : fusions)
        {
            final List<NeoEpitope> neDataList = Lists.newArrayList();

            // find all transcripts where the breakend is inside the coding region on the 5' gene
            final List<TranscriptData> upTransDataList = mGeneTransCache.getTranscripts(fusion.GeneIds[FS_UP]);
            final List<TranscriptData> downTransDataList = mGeneTransCache.getTranscripts(fusion.GeneIds[FS_DOWN]);
            final String[] validTranscriptNames = fusion.Transcripts;

            boolean sameGene = fusion.GeneIds[FS_UP].equals(fusion.GeneIds[FS_DOWN]);

            for(TranscriptData upTransData : upTransDataList)
            {
                if(!validTranscriptNames[FS_UP].contains(upTransData.TransName))
                    continue;

                if(!positionWithin(fusion.Positions[FS_UP], upTransData.CodingStart, upTransData.CodingEnd))
                    continue;

                for(TranscriptData downTransData : downTransDataList)
                {
                    // if the same gene then must be the same transcript
                    if(sameGene && upTransData.TransId != downTransData.TransId)
                        continue;

                    // must have a splice acceptor
                    if(downTransData.exons().size() <= 1)
                        continue;

                    int transRangeStart, transRangeEnd;

                    if(downTransData.Strand == POS_STRAND)
                    {
                        transRangeStart = downTransData.TransStart - DOWNSTREAM_PRE_GENE_DISTANCE;
                        transRangeEnd = downTransData.TransEnd;
                    }
                    else
                    {
                        transRangeStart = downTransData.TransStart;
                        transRangeEnd = downTransData.TransEnd + DOWNSTREAM_PRE_GENE_DISTANCE;
                    }

                    if(!positionWithin(fusion.Positions[FS_DOWN], transRangeStart, transRangeEnd))
                        continue;

                    if(!validTranscriptNames[FS_DOWN].contains(downTransData.TransName))
                        continue;

                    NeoEpitope neData = new SvNeoEpitope(fusion);
                    neDataList.add(neData);

                    neData.setTranscriptData(upTransData, downTransData);
                    neDataList.add(neData);
                }
            }

            processNeoEpitopes(neDataList);
        }
    }

    private void addPointMutations(final List<PointMutationData> pointMutations)
    {
        for(PointMutationData pointMutation : pointMutations)
        {
            final EnsemblGeneData geneData = mGeneTransCache.getGeneDataByName(pointMutation.Gene);

            final List<NeoEpitope> neDataList = Lists.newArrayList();

            final List<TranscriptData> transDataList = mGeneTransCache.getTranscripts(geneData.GeneId);

            for(TranscriptData transData : transDataList)
            {
                if(transData.CodingStart == null)
                    continue;

                if(!positionWithin(pointMutation.Position, transData.CodingStart, transData.CodingEnd))
                    continue;

                // check for a mutation within the stop codon at the bounds of the transcript
                if(pointMutation.Position <= transData.TransStart || pointMutation.Position >= transData.TransEnd)
                    continue;

                // must be exonic
                if(transData.exons().stream().noneMatch(x -> positionWithin(pointMutation.Position, x.Start, x.End)))
                    continue;

                if(transData.BioType.equals(BIOTYPE_NONSENSE_MED_DECAY))
                    continue;

                NeoEpitope neData = new PmNeoEpitope(pointMutation);
                neDataList.add(neData);

                neData.setTranscriptData(transData, transData);
            }

            processNeoEpitopes(neDataList);
        }
    }

    private void processNeoEpitopes(final List<NeoEpitope> neDataList)
    {
        if(neDataList.isEmpty())
            return;

        neDataList.forEach(x -> x.setCodingBases(mConfig.RefGenome, mConfig.RequiredAminoAcids));
        neDataList.forEach(x -> x.setAminoAcids());
        neDataList.forEach(x -> x.setNonsenseMediatedDecay());
        neDataList.forEach(x -> x.setSkippedSpliceSites(mGeneTransCache));

        // consolidate duplicates
        for(int i = 0; i < neDataList.size(); ++i)
        {
            NeoEpitope neData = neDataList.get(i);

            // filter out NEs with a novel stop codon
            if(neData.NovelAcid.equals(STOP_SYMBOL))
                continue;

            if(!neData.Valid)
            {
                IM_LOGGER.debug("skipping invalid neo: {}", neData);
                continue;
            }

            // filter out missense if in some transcripts their NE is the same as the wild-type
            if(neData.variantType() == NeoEpitopeType.MISSENSE && !neData.WildtypeAcids.isEmpty())
            {
                if(neData.aminoAcidString().contains(neData.WildtypeAcids))
                    continue;
            }

            final Set<String> upTransNames = Sets.newHashSet();
            final Set<String> downTransNames = Sets.newHashSet();
            upTransNames.add(neData.TransData[FS_UP].TransName);
            downTransNames.add(neData.TransData[FS_DOWN].TransName);

            final String aminoAcidStr = neData.aminoAcidString();

            if(!mConfig.WriteTransData)
            {
                int j = i + 1;
                while(j < neDataList.size())
                {
                    final NeoEpitope otherNeData = neDataList.get(j);
                    int minNmdCount = min(neData.NmdBasesMin, otherNeData.NmdBasesMin);
                    int maxNmdCount = max(neData.NmdBasesMax, otherNeData.NmdBasesMax);
                    int minScbCount = min(neData.StartCodonBasesMin, otherNeData.StartCodonBasesMin);
                    int maxScbCount = max(neData.StartCodonBasesMax, otherNeData.StartCodonBasesMax);

                    // remove exact matches or take the longer if one is a subset
                    if(aminoAcidStr.contains(otherNeData.aminoAcidString()))
                    {
                        neDataList.remove(j);
                    }
                    else if(otherNeData.aminoAcidString().contains(aminoAcidStr))
                    {
                        neDataList.set(i, otherNeData);
                        neData = otherNeData;
                        neDataList.remove(j);
                    }
                    else
                    {
                        ++j;
                        continue;
                    }

                    // take the shortest NMD base count
                    neData.NmdBasesMin = minNmdCount;
                    neData.NmdBasesMax = maxNmdCount;
                    neData.StartCodonBasesMin = minScbCount;
                    neData.StartCodonBasesMax = maxScbCount;

                    // collect up all transcripts
                    upTransNames.add(otherNeData.TransData[FS_UP].TransName);
                    downTransNames.add(otherNeData.TransData[FS_DOWN].TransName);
                }
            }

            writeData(neData, upTransNames, downTransNames);
        }
    }

    private void populateTpmMedians(
            final Set<String> upTransNames, final Set<String> downTransNames, final double[] tpmCancer, final double[] tpmCohort)
    {
        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            double cancerTotal = 0;
            double cohortTotal = 0;

            for(String transName : fs == FS_UP ? upTransNames : downTransNames)
            {
                final double[] result = mCohortTpmData.getTranscriptTpm(transName, mConfig.CancerType);
                cancerTotal += result[CANCER_VALUE];
                cohortTotal += result[COHORT_VALUE];
            }

            tpmCancer[fs] = cancerTotal;
            tpmCohort[fs] = cohortTotal;
        }
    }

    private void writeData(final NeoEpitope neData, final Set<String> upTransNames, final Set<String> downTransNames)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mWriter == null)
            {
                String outputFileName = mConfig.isMultiSample() ?
                        mConfig.OutputDir + "LNX_NEO_EPITOPES.csv" : mConfig.OutputDir + mCurrentSampleId + ".imu.neo_epitopes.csv";

                mWriter = createBufferedWriter(outputFileName, false);

                if(mConfig.isMultiSample())
                    mWriter.write("SampleId,");

                mWriter.write(NeoEpitopeFile.header());
                mWriter.newLine();
            }

            if(mConfig.isMultiSample())
                mWriter.write(String.format("%s,", mCurrentSampleId));

            final double[] tpmCancer = {0, 0};
            final double[] tpmCohort = {0, 0};
            populateTpmMedians(upTransNames, downTransNames, tpmCancer, tpmCohort);

            final NeoEpitopeFile neFile = neData.toFile(upTransNames, downTransNames, tpmCancer, tpmCohort);
            mWriter.write(NeoEpitopeFile.toString(neFile));
            mWriter.newLine();
        }
        catch (final IOException e)
        {
            IM_LOGGER.error("error writing neo-epitope output file: {}", e.toString());
        }
    }

    private void loadSvNeoEpitopes(final String filename)
    {
        if(filename == null || filename.isEmpty())
            return;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null)
            {
                IM_LOGGER.error("empty Linx neo-epitope file({})", filename);
                return;
            }

            boolean hasSampleId = line.contains(NE_SAMPLE_ID);
            final String singleSampleId = mConfig.SampleIds.size() == 1 ? mConfig.SampleIds.get(0) : "";

            int neCount = 0;
            String currentSampleId = "";
            List<NeoEpitopeFusion> fusions = null;

            while ((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(DELIMITER, -1);
                final String sampleId = hasSampleId ? items[0] : singleSampleId;

                if(!mConfig.SampleIds.contains(sampleId))
                    continue;

                if(!currentSampleId.equals(sampleId))
                {
                    fusions = Lists.newArrayList();
                    currentSampleId = sampleId;
                    mSampleFusionMap.put(sampleId, fusions);
                }

                fusions.add(NeoEpitopeFusion.fromString(line, hasSampleId));
                ++neCount;
            }

            IM_LOGGER.info("loaded {} Linx neo-epitope candidates from file: {}", neCount, filename);
        }
        catch(IOException exception)
        {
            IM_LOGGER.error("failed to read Linx neo-epitope file({})", filename, exception.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        NeoConfig.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        if(cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        NeoEpitopeAnnotator neoEpitopeAnnotator = new NeoEpitopeAnnotator(cmd);
        neoEpitopeAnnotator.run();

        IM_LOGGER.info("Neo-epitope annotations complete");
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
