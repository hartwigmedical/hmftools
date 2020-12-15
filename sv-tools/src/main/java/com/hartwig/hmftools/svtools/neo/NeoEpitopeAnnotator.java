package com.hartwig.hmftools.svtools.neo;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFusion.DELIMITER;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFusion.NE_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionWithin;
import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;
import static com.hartwig.hmftools.svtools.common.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.svtools.neo.NeoConfig.COHORT_FUSION_FILE;
import static com.hartwig.hmftools.svtools.neo.NeoConfig.GENE_TRANSCRIPTS_DIR;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.neo.NeoEpitopeFusion;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;
import org.jooq.Record;
import org.jooq.Result;

public class NeoEpitopeAnnotator
{
    public static final Logger IM_LOGGER = LogManager.getLogger(NeoEpitopeAnnotator.class);

    private final NeoConfig mConfig;
    private final Map<String,List<NeoEpitopeFusion>> mSampleFusionMap;

    private final EnsemblDataCache mGeneTransCache;
    private final DatabaseAccess mDbAccess;

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

        mDbAccess = createDatabaseAccess(cmd);

        loadSvNeoepitopes(cmd.getOptionValue(COHORT_FUSION_FILE));
    }

    public void run()
    {
        // check required inputs and config
        for(final String sampleId : mConfig.SampleIds)
        {
            final List<NeoEpitopeFusion> fusions = getSvFusions(sampleId);
            final List<PointMutationData> pointMutations = getSomaticVariants(sampleId);

            IM_LOGGER.debug("sample({}) loaded {} fusions and {} point mutations",
                    sampleId, fusions.size(), pointMutations.size());

            final List<NeoEpitope> neDataList = Lists.newArrayList();

            addSvFusions(fusions, neDataList);
            addPointMutations(pointMutations, neDataList);

            neDataList.forEach(x -> x.setCodingBases(mConfig.RefGenome, mConfig.RequiredAminoAcids));
            neDataList.forEach(x -> x.setAminoAcids());
            neDataList.forEach(x -> writeData(sampleId, x));
        }

        closeBufferedWriter(mWriter);
    }

    private final List<PointMutationData> getSomaticVariants(final String sampleId)
    {
        final List<PointMutationData> pointMutations = Lists.newArrayList();

        final Result<Record> result = mDbAccess.context().select().from(SOMATICVARIANT)
                .where(SOMATICVARIANT.SAMPLEID.eq(sampleId))
                .and(SOMATICVARIANT.FILTER.eq(PASS_FILTER))
                .and(SOMATICVARIANT.GENE.notEqual(""))
                .fetch();

        for (Record record : result)
        {
            final String gene = record.getValue(SOMATICVARIANT.GENE);

            if(gene.isEmpty() || mGeneTransCache.getGeneDataByName(gene) == null)
                continue;

            CodingEffect codingEffect = CodingEffect.valueOf(record.getValue(SOMATICVARIANT.WORSTCODINGEFFECT));

            if(codingEffect != NONSENSE_OR_FRAMESHIFT && codingEffect != MISSENSE)
                continue;

            String chromosome = record.getValue(SOMATICVARIANT.CHROMOSOME);
            int position = record.getValue(SOMATICVARIANT.POSITION);
            String ref = record.getValue(SOMATICVARIANT.REF);
            String alt = record.getValue(SOMATICVARIANT.ALT);
            double copyNumber = record.getValue(SOMATICVARIANT.VARIANTCOPYNUMBER);
            Integer localPhaseSet = record.get(SOMATICVARIANT.LOCALPHASESET);

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

    private void addSvFusions(final List<NeoEpitopeFusion> fusions, final List<NeoEpitope> neDataList)
    {
        for(NeoEpitopeFusion fusion : fusions)
        {
            // find all transcripts where the breakend is inside the coding region on the 5' gene
            final List<TranscriptData> upTransDataList = mGeneTransCache.getTranscripts(fusion.GeneIds[FS_UPSTREAM]);
            final List<TranscriptData> downTransDataList = mGeneTransCache.getTranscripts(fusion.GeneIds[FS_DOWNSTREAM]);

            for(TranscriptData upTransData : upTransDataList)
            {
                if(upTransData.CodingStart == null)
                    continue;

                if(!positionWithin(fusion.Positions[FS_UPSTREAM], upTransData.TransStart, upTransData.TransEnd))
                    continue;

                for(TranscriptData downTransData : downTransDataList)
                {
                    if(!positionWithin(fusion.Positions[FS_DOWNSTREAM], downTransData.TransStart, downTransData.TransEnd))
                        continue;

                    NeoEpitope neData = new SvNeoEpitope(fusion);
                    neDataList.add(neData);

                    neData.setTranscriptData(upTransData, downTransData);
                    neDataList.add(neData);
                }
            }
        }
    }

    private void addPointMutations(final List<PointMutationData> pointMutations, final List<NeoEpitope> neDataList)
    {
        for(PointMutationData pointMutation : pointMutations)
        {
            final EnsemblGeneData geneData = mGeneTransCache.getGeneDataByName(pointMutation.Gene);

            final List<TranscriptData> transDataList = mGeneTransCache.getTranscripts(geneData.GeneId);

            for(TranscriptData transData : transDataList)
            {
                if(transData.CodingStart == null)
                    continue;

                if(!positionWithin(pointMutation.Position, transData.CodingStart, transData.CodingEnd))
                    continue;

                NeoEpitope neData = new PmNeoEpitope(pointMutation);
                neDataList.add(neData);

                neData.setTranscriptData(transData, transData);
            }
        }
    }

    private boolean isDuplicate(final NeoEpitope neData)
    {
        return false;

        /*
        return mNeoEpitopeResults.stream()
                .anyMatch(x -> x.fusion().svId(true) == fusion.svId(true)
                        && x.fusion().svId(false) == fusion.svId(false));
        */
    }

    private void writeData(final String sampleId, final NeoEpitope neData)
    {
        /*
        if(mOutputDir.isEmpty())
            return;

        try
        {
            if(mFileWriter == null)
            {
                String outputFileName = mOutputDir + "LNX_NEO_EPITOPES.csv";

                mFileWriter = createBufferedWriter(outputFileName, false);

                mFileWriter.write("SampleId,Fusion,SameGene");
                mFileWriter.write(",UpstreamAminoAcids,DownstreamAminoAcids,NovelAminoAcid,NMDBases");

                for(int se = SE_START; se <= SE_END; ++se)
                {
                    String upDown = se == SE_START ? "Up" : "Down";

                    String fieldsStr = ",SvId" + upDown;
                    fieldsStr += ",Chr" + upDown;
                    fieldsStr += ",Pos" + upDown;
                    fieldsStr += ",Orient" + upDown;
                    fieldsStr += ",Trans" + upDown;
                    fieldsStr += ",Strand" + upDown;
                    fieldsStr += ",RegionType" + upDown;
                    fieldsStr += ",CodingType" + upDown;
                    fieldsStr += ",Exon" + upDown;
                    fieldsStr += ",Phase" + upDown;
                    mFileWriter.write(fieldsStr);
                }

                mFileWriter.newLine();
            }

            mFileWriter.write(String.format("%s,%s,%s",
                    sampleId, fusion.name(), fusion.upstreamTrans().geneName().equals(fusion.downstreamTrans().geneName())));

            mFileWriter.write(String.format(",%s,%s,%s,%d",
                    data.upstreamAcids(), data.downstreamAcids(), data.novelAcid(), data.downstreamNmdBases()));

            for(int se = SE_START; se <= SE_END; ++se)
            {
                boolean isUpstream = (se == SE_START);
                final Transcript trans = isUpstream ? fusion.upstreamTrans() : fusion.downstreamTrans();
                final GeneAnnotation gene = trans.gene();

                mFileWriter.write(String.format(",%d,%s,%d,%d",
                        gene.id(), gene.chromosome(), gene.position(), gene.orientation()));

                mFileWriter.write(String.format(",%s,%d,%s,%s,%d,%d",
                        trans.StableId, gene.Strand, trans.regionType(), trans.codingType(),
                        trans.nextSpliceExonRank(), trans.nextSpliceExonPhase()));
            }

            mFileWriter.newLine();
        }
        catch (final IOException e)
        {
            IM_LOGGER.error("error writing kataegis output file: {}", e.toString());
        }

        */
    }

    private void loadSvNeoepitopes(final String filename)
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

            int neCount = 0;
            String currentSampleId = "";
            List<NeoEpitopeFusion> fusions = null;

            while ((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(DELIMITER, -1);
                final String sampleId = hasSampleId ? items[0] : "";

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
