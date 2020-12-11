package com.hartwig.hmftools.svtools.neo;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFusion.DELIMITER;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFusion.NE_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.linx.fusion.FusionFinder.validFusionTranscript;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;
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
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.fusion.GeneAnnotation;
import com.hartwig.hmftools.common.fusion.Transcript;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.neo.NeoEpitopeFusion;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
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

            final List<NeoEpitopeData> neDataList = Lists.newArrayList();

            addSvFusions(fusions, neDataList);
            addPointMutations(pointMutations, neDataList);

            neDataList.forEach(x -> writeData(sampleId, x));
        }

        closeBufferedWriter(mWriter);
    }

    private final List<PointMutationData> getSomaticVariants(final String sampleId)
    {
        final List<PointMutationData> pointMutations = Lists.newArrayList();

        final List<SomaticVariant> somaticVariants = mDbAccess.readSomaticVariants(sampleId, VariantType.UNDEFINED);

        // filter to specific gene list
        for(final SomaticVariant variant : somaticVariants)
        {
            if(variant.gene().isEmpty() || mGeneTransCache.getGeneDataByName(variant.gene()) == null)
                continue;

            if(!variant.filter().equals(PASS_FILTER))
                continue;

            pointMutations.add(new PointMutationData(
                    variant.chromosome(), (int)variant.position(), variant.ref(), variant.alt(),
                    variant.gene(),
                    variant.localPhaseSet() != null ? variant.localPhaseSet() : -1));
        }

        return pointMutations;
    }

    private final List<NeoEpitopeFusion> getSvFusions(final String sampleId)
    {
        final List<NeoEpitopeFusion> fusions = mSampleFusionMap.get(sampleId);
        mSampleFusionMap.remove(sampleId);
        return fusions != null ? fusions : Lists.newArrayList();
    }

    private void addSvFusions(final List<NeoEpitopeFusion> fusions, final List<NeoEpitopeData> neDataList)
    {
        for(NeoEpitopeFusion fusion : fusions)
        {

        }

    }

    private void addPointMutations(final List<PointMutationData> pointMutations, final List<NeoEpitopeData> neDataList)
    {

    }

    // public final List<NeoEpitopeFusion> getResults() { return mNeoEpitopeResults; }

    /*
    public void reportNeoEpitopes(final String sampleId, final List<GeneFusion> fusions)
    {
        // mNeoEpitopeResults.clear();

        for(final GeneFusion fusion : fusions)
        {
            if(isDuplicate(fusion))
                continue;

            final Transcript upTrans = fusion.upstreamTrans();
            final Transcript downTrans = fusion.downstreamTrans();

            if(!upTrans.isCanonical() || !downTrans.isCanonical())
                continue;

            if(!validTranscripts(upTrans, downTrans))
                continue;

            boolean isPhased = fusion.phaseMatched()
                    && fusion.getExonsSkipped(true) == 0 && fusion.getExonsSkipped(false) == 0;

            int upstreamPhaseOffset = upTrans.ExonUpstreamPhase;
            int downstreamPhaseOffset = upstreamPhaseOffset == 0 || !isPhased ? 0 : 3 - upstreamPhaseOffset;

            IM_LOGGER.debug("fusion({}) SVs({} & {}) phased({}) upPhaseOffset({}) downPhaseOffset({})",
                    fusion.name(), fusion.svId(true), fusion.svId(false),
                    isPhased, upstreamPhaseOffset, downstreamPhaseOffset);


            String upstreamBases = getBaseString(
                    mConfig.RefGenome, upTrans, getTranscriptData(upTrans), AMINO_ACID_REF_COUNT, false, upstreamPhaseOffset);

            final TranscriptData downTransData = getTranscriptData(downTrans);

            String downstreamBases = getBaseString(
                    mConfig.RefGenome, downTrans, downTransData, AMINO_ACID_REF_COUNT, !isPhased, downstreamPhaseOffset);

            int nmdBaseCount = calcNonMediatedDecayBases(downTrans.gene(), downTransData);

            // upstream strand 1, bases will be retrieved from left to right (lower to higher), no need for any conversion
            // downstream strand 1, bases will be retrieved from left to right (lower to higher), no need for any conversion
            // upstream strand -1, bases will be retrieved from left to right (lower to higher), need to reverse and convert
            // downstream strand -1, bases will be retrieved from left to right (lower to higher), need to reverse and convert

            int upStrand = upTrans.gene().Strand;
            int downStrand = downTrans.gene().Strand;

            // correct for strand
            if(upStrand == -1)
                upstreamBases = reverseStrandBases(upstreamBases);

            if(downStrand == -1)
                downstreamBases = reverseStrandBases(downstreamBases);

            String novelCodonBases = "";

            if(upstreamPhaseOffset > upstreamBases.length() || downstreamPhaseOffset > downstreamBases.length())
                continue;

            // if upstream ends on a phase other than 0, need to take the bases from the downstream gene to make a novel codon

            if(upstreamPhaseOffset > 0)
            {
                // take the last 1 or 2 bases from the end of upstream gene's section
                novelCodonBases = upstreamBases.substring(upstreamBases.length() - upstreamPhaseOffset);
                upstreamBases = upstreamBases.substring(0, upstreamBases.length() - upstreamPhaseOffset);
            }

            if(isPhased)
            {
                novelCodonBases += downstreamBases.substring(0, downstreamPhaseOffset);
                downstreamBases = downstreamBases.substring(downstreamPhaseOffset);
            }
            else
            {
                novelCodonBases += downstreamBases;
                downstreamBases = "";
            }

            IM_LOGGER.debug("ne({}) upBases({}) novelCodon({}) downBases({}) downNmdBases({})",
                    neData, upstreamBases, checkTrimBases(novelCodonBases), checkTrimBases(downstreamBases), nmdBaseCount);

            if(upstreamBases.isEmpty() || (downstreamBases.isEmpty() && novelCodonBases.isEmpty()))
                continue;

            final String upstreamRefAminoAcids = getAminoAcids(upstreamBases, false);
            final String novelAminoAcids = getAminoAcids(novelCodonBases, !isPhased);
            String downstreamRefAminoAcids = getAminoAcids(downstreamBases, !isPhased);

            IM_LOGGER.debug("fusion({}) upAA({}) novel({}) downAA({})",
                    fusion.name(), upstreamRefAminoAcids, checkTrimBases(novelAminoAcids), checkTrimBases(downstreamRefAminoAcids));

            if(novelAminoAcids.equals(STOP_SYMBOL))
                downstreamRefAminoAcids = "";

            NeoEpitopeFusion neoEpData = ImmutableNeoEpitopeData.builder()
                    .fusion(fusion)
                    .upstreamAcids(upstreamRefAminoAcids)
                    .downstreamAcids(downstreamRefAminoAcids)
                    .novelAcid(novelAminoAcids)
                    .downstreamNmdBases(nmdBaseCount)
                    .build();

            mNeoEpitopeResults.add(neoEpData);
        }
    }
    */

    private boolean isDuplicate(final NeoEpitopeData neData)
    {
        return false;

        /*
        return mNeoEpitopeResults.stream()
                .anyMatch(x -> x.fusion().svId(true) == fusion.svId(true)
                        && x.fusion().svId(false) == fusion.svId(false));
        */
    }



    private TranscriptData getTranscriptData(final Transcript transcript)
    {
        final TranscriptData transData = mGeneTransCache.getTranscriptData(transcript.gene().StableId, transcript.StableId);

        if(transData == null)
        {
            IM_LOGGER.error("gene({}) transcript({}) data not found", transcript.gene().GeneName, transcript.StableId);
            return null;
        }

        return transData;
    }

    /*

    public void checkFusions(List<GeneFusion> existingFusions, List<GeneAnnotation> genesList1,List<GeneAnnotation> genesList2)
    {
        for (final GeneAnnotation gene1 : genesList1)
        {
            boolean startUpstream = gene1.isUpstream();

            final Transcript trans1 = gene1.transcripts().stream().filter(Transcript::isCanonical).findFirst().orElse(null);

            if (trans1 == null)
                continue;

            for (final GeneAnnotation gene2 : genesList2)
            {
                boolean endUpstream = gene2.isUpstream();

                if (startUpstream == endUpstream)
                    continue;

                final Transcript trans2 = gene2.transcripts().stream().filter(Transcript::isCanonical).findFirst().orElse(null);

                if (trans2 == null)
                    continue;

                final Transcript upstreamTrans = startUpstream ? trans1 : trans2;
                final Transcript downstreamTrans = upstreamTrans == trans1 ? trans2 : trans1;

                if(existingFusions.stream().anyMatch(x -> x.upstreamTrans() == upstreamTrans && x.downstreamTrans() == downstreamTrans))
                    continue;

                GeneFusion newFusion = checkNeoEpitopeFusion(upstreamTrans, downstreamTrans);

                if(newFusion != null)
                    existingFusions.add(newFusion);
            }
        }
    }

     */

    private boolean validTranscripts(final Transcript upstreamTrans, final Transcript downstreamTrans)
    {
        if(!validFusionTranscript(upstreamTrans))
            return false;

        if(downstreamTrans.nonCoding() || downstreamTrans.ExonMax == 1)
            return false;

        if(upstreamTrans.preCoding() || upstreamTrans.ExonUpstreamPhase == -1)
            return false;

        if(upstreamTrans.isExonic() && downstreamTrans.isExonic()) // rare and too complicated for now
            return false;

        return true;
    }

    private void writeData(final String sampleId, final NeoEpitopeData neData)
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
