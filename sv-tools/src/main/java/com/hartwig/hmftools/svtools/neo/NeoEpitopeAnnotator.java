package com.hartwig.hmftools.svtools.neo;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.switchStream;
import static com.hartwig.hmftools.common.fusion.Transcript.CODING_BASES;
import static com.hartwig.hmftools.common.fusion.Transcript.TOTAL_CODING_BASES;
import static com.hartwig.hmftools.common.fusion.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.fusion.TranscriptCodingType.NON_CODING;
import static com.hartwig.hmftools.common.fusion.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.fusion.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.fusion.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.fusion.TranscriptRegionType.INTRONIC;
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
import static com.hartwig.hmftools.svtools.neo.AminoAcidConverter.STOP_SYMBOL;
import static com.hartwig.hmftools.svtools.neo.AminoAcidConverter.reverseStrandBases;
import static com.hartwig.hmftools.svtools.neo.NeoConfig.AMINO_ACID_REF_COUNT;
import static com.hartwig.hmftools.svtools.neo.NeoConfig.COHORT_FUSION_FILE;
import static com.hartwig.hmftools.svtools.neo.NeoConfig.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.svtools.neo.NeoUtils.calcNonMediatedDecayBases;
import static com.hartwig.hmftools.svtools.neo.NeoUtils.checkTrimBases;
import static com.hartwig.hmftools.svtools.neo.NeoUtils.getAminoAcids;
import static com.hartwig.hmftools.svtools.neo.NeoUtils.getCodingBases;

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
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.fusion.Transcript;
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

            final List<NeoEpitopeData> neDataList = Lists.newArrayList();

            addSvFusions(fusions, neDataList);
            addPointMutations(pointMutations, neDataList);

            neDataList.forEach(x -> setContextBases(x));
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

    private void addSvFusions(final List<NeoEpitopeFusion> fusions, final List<NeoEpitopeData> neDataList)
    {
        for(NeoEpitopeFusion fusion : fusions)
        {
            // find all transcripts where the breakend is inside the coding region on the 5' gene
            for(int fs = FS_UPSTREAM; fs <= FS_DOWNSTREAM; ++fs)
            {
                final List<TranscriptData> transDataList = mGeneTransCache.getTranscripts(fusion.GeneIds[fs]);

                for(TranscriptData transData : transDataList)
                {
                    if(fs == FS_UPSTREAM && transData.CodingStart == null)
                        continue;

                    if(!positionWithin(fusion.Positions[fs], transData.TransStart, transData.TransEnd))
                        continue;

                    NeoEpitopeData neData = new NeoEpitopeData(null, fusion);

                    neData.TransData[fs] = transData;

                    setTranscriptContext(neData, transData, fusion.Positions[fs], fs);
                    setTranscriptCodingData(neData, transData, fusion.Positions[fs], fusion.InsertSequence.length(), fs);
                }
            }
        }
    }

    private void addPointMutations(final List<PointMutationData> pointMutations, final List<NeoEpitopeData> neDataList)
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

                NeoEpitopeData neData = new NeoEpitopeData(pointMutation, null);

                neData.TransData[FS_UPSTREAM] = neData.TransData[FS_DOWNSTREAM] = transData;

                int indelBaseDiff = pointMutation.Alt.length() - pointMutation.Ref.length();
                int insertedBaseLength = max(indelBaseDiff, 0);

                // set the data for the lower part of the mutation
                int lowerStream = geneData.Strand == POS_STRAND ? FS_UPSTREAM : FS_DOWNSTREAM;

                setTranscriptContext(neData, transData, pointMutation.Position, lowerStream);

                int insSeqLength = lowerStream == FS_UPSTREAM ? insertedBaseLength : 0;

                setTranscriptCodingData(neData, transData, pointMutation.Position, insSeqLength, lowerStream);

                int upperStream = switchStream(lowerStream);

                // for DELs, set the downstream data as well since it can cross exon-boundaries and/or affect coding bases
                if(indelBaseDiff != 0)
                {
                    int adjustedPosition = indelBaseDiff < 0 ? pointMutation.Position + abs(indelBaseDiff) : pointMutation.Position;
                    setTranscriptContext(neData, transData, adjustedPosition, upperStream);

                    insSeqLength = upperStream == FS_UPSTREAM ? insertedBaseLength : 0;
                    setTranscriptCodingData(neData, transData, adjustedPosition, insSeqLength, upperStream);
                }
                else
                {
                    neData.RegionType[upperStream] = neData.RegionType[lowerStream];
                    neData.CodingType[upperStream] = neData.CodingType[lowerStream];
                    neData.Phases[upperStream] = neData.Phases[lowerStream];

                    if(neData.RegionType[lowerStream] == INTRONIC)
                    {
                        if(lowerStream == FS_UPSTREAM)
                            neData.ExonRank[upperStream] = neData.ExonRank[lowerStream] + 1;
                        else
                            neData.ExonRank[upperStream] = neData.ExonRank[lowerStream] - 1;
                    }
                }
            }
        }
    }

    private void setTranscriptContext(
            final NeoEpitopeData neData, final TranscriptData transData, int position, int fs)
    {
        // determine phasing, coding and region context
        boolean isUpstream = fs == FS_UPSTREAM;

        for(ExonData exon : transData.exons())
        {
            if(position > exon.ExonEnd)
                continue;

            if(position < exon.ExonStart)
            {
                // intronic
                neData.RegionType[fs] = INTRONIC;

                // upstream pos-strand, before next exon then take the exon before's rank
                if(transData.Strand == POS_STRAND && isUpstream)
                    neData.ExonRank[fs] = exon.ExonRank - 1;
                else if(transData.Strand == NEG_STRAND && isUpstream)
                    neData.ExonRank[fs] = exon.ExonRank;
                else if(transData.Strand == POS_STRAND && !isUpstream)
                    neData.ExonRank[fs] = exon.ExonRank;
                else if(transData.Strand == NEG_STRAND && !isUpstream)
                    neData.ExonRank[fs] = exon.ExonRank - 1;

                neData.Phases[fs] = exon.ExonPhase;
            }
            else if(positionWithin(position, exon.ExonStart, exon.ExonEnd))
            {
                neData.RegionType[fs] = EXONIC;
                neData.ExonRank[fs] = exon.ExonRank;
            }
        }
    }

    private void setTranscriptCodingData(
            final NeoEpitopeData neData, final TranscriptData transData, int position, int insSeqLength, int fs)
    {
        if(transData.CodingStart != null)
        {
            if(positionWithin(position, transData.CodingStart, transData.CodingEnd))
                neData.CodingType[fs] = CODING;
            else if(transData.Strand == POS_STRAND && position < transData.CodingStart)
                neData.CodingType[fs] = UTR_5P;
            else if(transData.Strand == NEG_STRAND && position > transData.CodingEnd)
                neData.CodingType[fs] = UTR_5P;
            else
                neData.CodingType[fs] = UTR_3P;
        }
        else
        {
            neData.CodingType[fs] = NON_CODING;
            neData.Phases[fs] = -1;
        }

        if(neData.CodingType[fs] == CODING && neData.RegionType[fs] == EXONIC)
        {
            int[] codingData = Transcript.calcCodingBases(transData.CodingStart, transData.CodingEnd, transData.exons(), position);
            int codingBases = transData.Strand == POS_STRAND ? codingData[CODING_BASES] : codingData[TOTAL_CODING_BASES] - codingData[CODING_BASES];

            codingBases -= 1;

            // factor in insert sequence for the upstream partner
            codingBases += insSeqLength;

            neData.Phases[fs] = codingBases % 3;
        }
    }

    // public final List<NeoEpitopeFusion> getResults() { return mNeoEpitopeResults; }

    private void setContextBases(final NeoEpitopeData neData)
    {
        boolean isPhased = neData.phaseMatched();

        int upstreamPhaseOffset = neData.Phases[FS_UPSTREAM];
        int downstreamPhaseOffset = upstreamPhaseOffset == 0 || !isPhased ? 0 : 3 - upstreamPhaseOffset;

        IM_LOGGER.debug("fusion({}) phased({}) upPhaseOffset({}) downPhaseOffset({})",
                neData.toString(), isPhased, upstreamPhaseOffset, downstreamPhaseOffset);

        String upstreamBases = getCodingBases(mConfig.RefGenome, neData, FS_UPSTREAM, AMINO_ACID_REF_COUNT, upstreamPhaseOffset);

        String downstreamBases = getCodingBases(mConfig.RefGenome, neData, FS_DOWNSTREAM, AMINO_ACID_REF_COUNT, downstreamPhaseOffset);

        neData.DownstreamNmdBases = calcNonMediatedDecayBases(neData, FS_DOWNSTREAM);

        // upstream strand 1, bases will be retrieved from left to right (lower to higher), no need for any conversion
        // downstream strand 1, bases will be retrieved from left to right (lower to higher), no need for any conversion
        // upstream strand -1, bases will be retrieved from left to right (lower to higher), need to reverse and convert
        // downstream strand -1, bases will be retrieved from left to right (lower to higher), need to reverse and convert

        // correct for strand
        if(neData.strand(FS_UPSTREAM) == NEG_STRAND)
            upstreamBases = reverseStrandBases(upstreamBases);

        if(neData.strand(FS_DOWNSTREAM) == NEG_STRAND)
            downstreamBases = reverseStrandBases(downstreamBases);

        String novelCodonBases = "";

        if(upstreamPhaseOffset > upstreamBases.length() || downstreamPhaseOffset > downstreamBases.length())
        {
            IM_LOGGER.error("ne({}) invalid upBases({} phaseOffset={}) or downBases({} phaseOffset={})",
                    neData, upstreamBases, upstreamPhaseOffset, downstreamBases, downstreamPhaseOffset);
            return;
        }

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
                neData, upstreamBases, checkTrimBases(novelCodonBases), checkTrimBases(downstreamBases), neData.DownstreamNmdBases);

        final String upstreamRefAminoAcids = getAminoAcids(upstreamBases, false);
        final String novelAminoAcids = getAminoAcids(novelCodonBases, !isPhased);
        String downstreamRefAminoAcids = getAminoAcids(downstreamBases, !isPhased);

        IM_LOGGER.debug("ne({}) upAA({}) novel({}) downAA({})",
                neData, upstreamRefAminoAcids, checkTrimBases(novelAminoAcids), checkTrimBases(downstreamRefAminoAcids));

        if(novelAminoAcids.equals(STOP_SYMBOL))
            downstreamRefAminoAcids = "";

        neData.UpstreamAcids = upstreamRefAminoAcids;
        neData.DownstreamAcids = downstreamRefAminoAcids;
        neData.NovelAcid = novelAminoAcids;

    }

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
