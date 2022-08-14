package com.hartwig.hmftools.neo.epitope;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.AminoAcidRna.STOP_SYMBOL;
import static com.hartwig.hmftools.common.gene.TranscriptProteinData.BIOTYPE_NONSENSE_MED_DECAY;
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
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.localPhaseSetsStringToList;
import static com.hartwig.hmftools.neo.NeoCommon.DOWNSTREAM_PRE_GENE_DISTANCE;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.epitope.NeoEpitopeAnnotator.initialiseNeoepitopeWriter;
import static com.hartwig.hmftools.neo.epitope.NeoEpitopeAnnotator.initialisePeptideWriter;
import static com.hartwig.hmftools.neo.epitope.NeoEpitopeAnnotator.writeNeoepitopes;
import static com.hartwig.hmftools.neo.epitope.NeoEpitopeAnnotator.writePeptideHlaData;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Somaticvariant.SOMATICVARIANT;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
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
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.neo.NeoEpitopeFusion;
import com.hartwig.hmftools.common.neo.NeoEpitopeType;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;

import org.jooq.Record;
import org.jooq.Record8;
import org.jooq.Result;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.vcf.VCFCodec;

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
            final CohortTpmData cohortTpmData)
    {
        mSampleData = sampleData;

        mConfig = config;
        mGeneTransCache = ensemblDataCache;
        mDbAccess = dbAccess;
        mCohortTpmData = cohortTpmData;

        mNextNeoEpitopeId = 0;
        mFusions = Lists.newArrayList();

        mNeoEpitopeWriter = initialiseNeoepitopeWriter(mConfig.OutputDir, mSampleData.Id);
        mPeptideWriter = mConfig.WritePeptides ? initialisePeptideWriter(mConfig.OutputDir, mSampleData.Id) : null;
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

        if(mDbAccess != null)
        {
            final Result<Record8<String, String, String, Integer, String, String, Double, String>> result = mDbAccess.context()
                    .select(SOMATICVARIANT.GENE, SOMATICVARIANT.WORSTCODINGEFFECT, SOMATICVARIANT.CHROMOSOME, SOMATICVARIANT.POSITION,
                            SOMATICVARIANT.REF, SOMATICVARIANT.ALT, SOMATICVARIANT.VARIANTCOPYNUMBER, SOMATICVARIANT.LOCALPHASESET)
                    .from(Tables.SOMATICVARIANT)
                    .where(Tables.SOMATICVARIANT.SAMPLEID.eq(mSampleData.Id))
                    .and(Tables.SOMATICVARIANT.FILTER.eq(PASS_FILTER))
                    .and(Tables.SOMATICVARIANT.GENE.notEqual(""))
                    .and(Tables.SOMATICVARIANT.WORSTCODINGEFFECT.in(validCodingEffects))
                    .fetch();

            for(Record record : result)
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
                String localPhaseSetStr = record.get(Tables.SOMATICVARIANT.LOCALPHASESET);
                List<Integer> localPhaseSets = localPhaseSetsStringToList(localPhaseSetStr);

                pointMutations.add(new PointMutationData(chromosome, position, ref, alt, gene,
                        codingEffect, copyNumber, localPhaseSets != null ? localPhaseSets.get(0) : -1));
            }
        }
        else
        {
            if(mConfig.SomaticVcf == null)
                return pointMutations;

            String somaticVcf = mConfig.SomaticVcf.contains("*") ?
                    mConfig.SomaticVcf.replaceAll("\\*", mSampleData.Id) : mConfig.SomaticVcf;

            if(!Files.exists(Paths.get(somaticVcf)))
            {
                NE_LOGGER.warn("Purple somatic VCF file({}) not found", somaticVcf);
                return pointMutations;
            }

            try
            {
                CompoundFilter filter = new CompoundFilter(true);
                filter.add(new PassingVariantFilter());

                SomaticVariantFactory variantFactory = new SomaticVariantFactory(filter);

                final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(somaticVcf, new VCFCodec(), false);

                for (VariantContext variant : reader.iterator())
                {
                    if (filter.test(variant))
                    {
                        final SomaticVariant somaticVariant = variantFactory.createVariant(mSampleData.Id, variant).orElse(null);

                        if(somaticVariant == null)
                            continue;

                        if(somaticVariant.gene().isEmpty() || mGeneTransCache.getGeneDataByName(somaticVariant.gene()) == null)
                            continue;

                        if(somaticVariant.worstCodingEffect() != NONSENSE_OR_FRAMESHIFT && somaticVariant.worstCodingEffect() != MISSENSE)
                            continue;

                        pointMutations.add(new PointMutationData(
                                somaticVariant.chromosome(), (int)somaticVariant.position(), somaticVariant.ref(), somaticVariant.alt(),
                                somaticVariant.gene(), somaticVariant.worstCodingEffect(), somaticVariant.adjustedCopyNumber(),
                                somaticVariant.localPhaseSets() != null ? somaticVariant.topLocalPhaseSet() : -1));
                    }
                }

                NE_LOGGER.debug("loaded {} somatic variants from file({})", pointMutations.size(), somaticVcf);
            }
            catch(IOException e)
            {
                NE_LOGGER.error(" failed to read somatic VCF file({}): {}", somaticVcf, e.toString());
            }
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
            final GeneData geneData = mGeneTransCache.getGeneDataByName(pointMutation.Gene);

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
}
