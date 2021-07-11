package com.hartwig.hmftools.neo;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.AminoAcidConverter.STOP_SYMBOL;
import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_NONSENSE_MED_DECAY;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFusion.generateFilename;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.neo.NeoCommon.DOWNSTREAM_PRE_GENE_DISTANCE;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.NeoEpitopeAnnotator.initialiseNeoepitopeWriter;
import static com.hartwig.hmftools.neo.NeoEpitopeAnnotator.initialisePeptideWriter;
import static com.hartwig.hmftools.neo.NeoEpitopeAnnotator.writeNeoepitopes;
import static com.hartwig.hmftools.neo.NeoEpitopeAnnotator.writePeptideHlaData;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.neo.NeoEpitopeFusion;
import com.hartwig.hmftools.common.neo.NeoEpitopeType;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;

import org.jooq.Record;
import org.jooq.Result;

public class NeoSampleTask implements Callable
{
    private final SampleData mSampleData;

    private final NeoConfig mConfig;
    private final List<NeoEpitopeFusion> mFusions;

    private final EnsemblDataCache mGeneTransCache;
    private final DatabaseAccess mDbAccess;
    private final CohortTpmData mCohortTpmData;

    private int mNextNeoEpitopeId;
    private final BufferedWriter mNeoEpitopeWriter;
    private BufferedWriter mPeptideWriter;

    public NeoSampleTask(
            final SampleData sampleData, final NeoConfig config, final EnsemblDataCache ensemblDataCache, final DatabaseAccess dbAccess,
            final CohortTpmData cohortTpmData, final BufferedWriter neoEpitopeWriter, final BufferedWriter peptideWriter)
    {
        mSampleData = sampleData;

        mConfig = config;
        mGeneTransCache = ensemblDataCache;
        mDbAccess = dbAccess;
        mCohortTpmData = cohortTpmData;

        mNextNeoEpitopeId = 0;
        mFusions = Lists.newArrayList();

        if(neoEpitopeWriter != null && peptideWriter != null)
        {
            mNeoEpitopeWriter = null;
            mPeptideWriter = null;
        }
        else
        {
            mNeoEpitopeWriter = initialiseNeoepitopeWriter(mConfig.OutputDir, mSampleData.Id);
            mPeptideWriter = initialisePeptideWriter(mConfig.OutputDir, mSampleData.Id);
        }
    }

    @Override
    public Long call()
    {
        processSample();
        return (long)1;
    }

    public void processSample()
    {
        final List<NeoEpitopeFusion> fusions = getSvFusions();
        final List<PointMutationData> pointMutations = getSomaticVariants();

        NE_LOGGER.info("sample({}) loaded {} fusions and {} point mutations",
                mSampleData.Id, fusions.size(), pointMutations.size());

        addSvFusions(fusions);
        addPointMutations(pointMutations);

        closeBufferedWriter(mNeoEpitopeWriter);
        closeBufferedWriter(mPeptideWriter);
    }

    private final List<PointMutationData> getSomaticVariants()
    {
        List<PointMutationData> pointMutations = Lists.newArrayList();

        List<String> validCodingEffects = Lists.newArrayList(NONSENSE_OR_FRAMESHIFT.toString(), MISSENSE.toString());

        final Result<Record> result = mDbAccess.context().select().from(Tables.SOMATICVARIANT)
                .where(Tables.SOMATICVARIANT.SAMPLEID.eq(mSampleData.Id))
                .and(Tables.SOMATICVARIANT.FILTER.eq(PASS_FILTER))
                .and(Tables.SOMATICVARIANT.GENE.notEqual(""))
                .and(Tables.SOMATICVARIANT.WORSTCODINGEFFECT.in(validCodingEffects))
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

    private List<NeoEpitopeFusion> getSvFusions()
    {
        List<NeoEpitopeFusion> fusions = Lists.newArrayList();

        if(mConfig.SvFusionsDir == null)
            return fusions;

        final String filename = generateFilename(mConfig.SvFusionsDir, mSampleData.Id);

        if(!Files.exists(Paths.get(filename)))
        {
            NE_LOGGER.warn("Linx neo-epitope file({}) not found", filename);
            return fusions;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null)
            {
                NE_LOGGER.error("empty Linx neo-epitope file({})", filename);
                return fusions;
            }

            int neCount = 0;

            while ((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(DELIMITER, -1);

                fusions.add(NeoEpitopeFusion.fromString(line, false));
                ++neCount;
            }

            NE_LOGGER.debug("loaded {} Linx neo-epitope candidates from file: {}", neCount, filename);
        }
        catch(IOException exception)
        {
            NE_LOGGER.error("failed to read Linx neo-epitope file({})", filename, exception.toString());
        }

        return fusions;
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
        neDataList.forEach(x -> x.setAminoAcids(mConfig.RefGenome, mConfig.RequiredAminoAcids));
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
                NE_LOGGER.debug("skipping invalid neo: {}", neData);
                continue;
            }

            // filter out missense if in some transcripts their NE is the same as the wild-type
            if(neData.variantType() == NeoEpitopeType.MISSENSE && !neData.wildtypeAcids().isEmpty())
            {
                if(neData.aminoAcidString().contains(neData.wildtypeAcids()))
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
                    int minCodingBaseLen = min(neData.CodingBasesLengthMin, otherNeData.CodingBasesLengthMin);
                    int maxCodingBaseLen = max(neData.CodingBasesLengthMax, otherNeData.CodingBasesLengthMax);

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
                    neData.CodingBasesLengthMin = minCodingBaseLen;
                    neData.CodingBasesLengthMax = maxCodingBaseLen;

                    // collect up all transcripts
                    upTransNames.add(otherNeData.TransData[FS_UP].TransName);
                    downTransNames.add(otherNeData.TransData[FS_DOWN].TransName);
                }
            }

            int neId = mNextNeoEpitopeId++;
            writeNeoepitopes(mNeoEpitopeWriter, mSampleData, false, mCohortTpmData, neId, neData, upTransNames, downTransNames);
            writePeptideHlaData(mPeptideWriter, mSampleData, false, mConfig, neId, neData);
        }
    }

    /*
    private void writeData(int neId, final NeoEpitope neData, final Set<String> upTransNames, final Set<String> downTransNames)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mConfig.WriteCohortFile)
                mNeoEpitopeWriter.write(String.format("%s,", mSampleData.Id));

            final double[] tpmCancer = {0, 0};
            final double[] tpmCohort = {0, 0};
            populateTpmMedians(upTransNames, downTransNames, tpmCancer, tpmCohort);

            final NeoEpitopeFile neFile = neData.toFile(neId, upTransNames, downTransNames, tpmCancer, tpmCohort);
            mNeoEpitopeWriter.write(NeoEpitopeFile.toString(neFile));
            mNeoEpitopeWriter.newLine();
        }
        catch (final IOException e)
        {
            IM_LOGGER.error("error writing neo-epitope output file: {}", e.toString());
        }
    }

    private void writePeptideHlaData(int neId, final NeoEpitope neData)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        if(!mConfig.CommonHlaTypes.isEmpty() && mSampleData.HlaTypes.isEmpty())
        {
            mSampleData.HlaTypes.addAll(mConfig.CommonHlaTypes);
        }

        if(mSampleData.HlaTypes.isEmpty() || mConfig.PeptideLengths[SE_START] == 0 || mConfig.PeptideLengths[SE_END] == 0)
            return;

        try
        {
            if(mConfig.WriteCohortFile)
                mPeptideWriter.write(String.format("%s,", mSampleData.Id));

            final List<PeptideData> peptides = NeoUtils.generatePeptides(
                    neData.UpstreamAcids, neData.NovelAcid, neData.DownstreamAcids, mConfig.PeptideLengths, mConfig.PeptideFlanks);

            Set<String> uniqueHlaTypes = Sets.newHashSet();

            for(String hlaType : mSampleData.HlaTypes)
            {
                if(uniqueHlaTypes.contains(hlaType))
                    continue;

                uniqueHlaTypes.add(hlaType);

                final String predictionHlaType = NeoUtils.convertHlaTypeForPredictions(hlaType);

                if(predictionHlaType == null)
                {
                    IM_LOGGER.error("sample({} skipping invalid HLA type: {}", mSampleData.Id, hlaType);
                    continue;
                }

                for(PeptideData peptideData : peptides)
                {
                    // skip any peptide which is contained within the upstream wildtype AAs
                    if(neData.UpstreamWildTypeAcids.contains(peptideData.Peptide))
                        continue;

                    // for now skip any upstream peptide containing the 21st AA until MhcFlurry can handle it
                    if(peptideData.Peptide.contains(AA_SELENOCYSTEINE) || peptideData.UpFlank.contains(AA_SELENOCYSTEINE))
                        continue;

                    mPeptideWriter.write(String.format("%d,%s,%s,%s,%s",
                            neId, predictionHlaType, peptideData.Peptide, peptideData.UpFlank, peptideData.DownFlank));
                    mPeptideWriter.newLine();
                }
            }
        }
        catch (final IOException e)
        {
            IM_LOGGER.error("error writing HLA peptide output file: {}", e.toString());
        }
    }
    */

}
