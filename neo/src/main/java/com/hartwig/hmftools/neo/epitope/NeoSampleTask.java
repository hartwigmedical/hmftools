package com.hartwig.hmftools.neo.epitope;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.AminoAcidRna.STOP_SYMBOL;
import static com.hartwig.hmftools.common.gene.TranscriptProteinData.BIOTYPE_NONSENSE_MED_DECAY;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.ITEM_DELIM;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFusion.generateFilename;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.neo.NeoCommon.DOWNSTREAM_PRE_GENE_DISTANCE;
import static com.hartwig.hmftools.neo.NeoCommon.IMMUNE_TRANSCRIPT_PREFIXES;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.epitope.NeoEpitopeFinder.initialiseNeoepitopeWriter;
import static com.hartwig.hmftools.neo.epitope.NeoEpitopeFinder.writeNeoepitopes;
import static com.hartwig.hmftools.neo.epitope.PointMutationData.isRelevantMutation;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.BufferedWriter;
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
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.vcf.VCFCodec;

public class NeoSampleTask implements Callable
{
    private final String mSampleId;

    private final NeoConfig mConfig;

    private final EnsemblDataCache mGeneTransCache;

    private int mNextNeoEpitopeId;
    private final BufferedWriter mWriter;

    public NeoSampleTask(final String sampleId, final NeoConfig config, final EnsemblDataCache ensemblDataCache)
    {
        mSampleId = sampleId;

        mConfig = config;
        mGeneTransCache = ensemblDataCache;

        mNextNeoEpitopeId = 0;

        mWriter = initialiseNeoepitopeWriter(mConfig.OutputDir, mSampleId);
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
                mSampleId, fusions.size(), pointMutations.size());

        addSvFusions(fusions);
        addPointMutations(pointMutations);

        closeBufferedWriter(mWriter);
    }

    private List<PointMutationData> getSomaticVariants()
    {
        List<PointMutationData> pointMutations = Lists.newArrayList();

        if(mConfig.SomaticVcf == null)
            return pointMutations;

        String somaticVcf = mConfig.SomaticVcf.contains("*") ?
                mConfig.SomaticVcf.replaceAll("\\*", mSampleId) : mConfig.SomaticVcf;

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

            for(VariantContext variant : reader.iterator())
            {
                if(filter.test(variant))
                {
                    final SomaticVariant somaticVariant = variantFactory.createVariant(mSampleId, variant).orElse(null);

                    if(somaticVariant == null)
                        continue;

                    if(somaticVariant.gene().isEmpty() || mGeneTransCache.getGeneDataByName(somaticVariant.gene()) == null)
                        continue;

                    if(!isRelevantMutation(somaticVariant))
                        continue;

                    pointMutations.add(new PointMutationData(
                            somaticVariant.chromosome(), somaticVariant.position(), somaticVariant.ref(), somaticVariant.alt(),
                            somaticVariant.gene(), somaticVariant.worstCodingEffect(), somaticVariant.adjustedCopyNumber(),
                            somaticVariant.subclonalLikelihood(),
                            somaticVariant.localPhaseSets() != null ? somaticVariant.topLocalPhaseSet() : -1));
                }
            }

            NE_LOGGER.debug("loaded {} somatic variants from file({})", pointMutations.size(), somaticVcf);
        }
        catch(IOException e)
        {
            NE_LOGGER.error(" failed to read somatic VCF file({}): {}", somaticVcf, e.toString());
        }

        return pointMutations;
    }

    private List<NeoEpitopeFusion> getSvFusions()
    {
        List<NeoEpitopeFusion> fusions = Lists.newArrayList();

        if(mConfig.SvFusionsDir == null)
            return fusions;

        final String filename = generateFilename(mConfig.SvFusionsDir, mSampleId);

        if(!Files.exists(Paths.get(filename)))
        {
            NE_LOGGER.warn("Linx neo-epitope file({}) not found", filename);
            return fusions;
        }

        try
        {
            fusions.addAll(NeoEpitopeFusion.read(filename));
            NE_LOGGER.debug("loaded {} Linx neo-epitope candidates from file: {}", fusions.size(), filename);
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

            String[] upTransNames = fusion.Transcripts[FS_UP].split(ITEM_DELIM, -1);
            String[] downTransNames = fusion.Transcripts[FS_DOWN].split(ITEM_DELIM, -1);

            boolean sameGene = fusion.GeneIds[FS_UP].equals(fusion.GeneIds[FS_DOWN]);

            for(String upTransName : upTransNames)
            {
                TranscriptData upTransData = upTransDataList.stream().filter(x -> x.TransName.equals(upTransName)).findFirst().orElse(null);

                if(upTransData == null)
                    continue;

                if(!positionWithin(fusion.Positions[FS_UP], upTransData.CodingStart, upTransData.CodingEnd))
                    continue;

                for(String downTransName : downTransNames)
                {
                    TranscriptData downTransData = downTransDataList.stream()
                            .filter(x -> x.TransName.equals(downTransName)).findFirst().orElse(null);

                    if(downTransData == null)
                        continue;

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

                    NeoEpitope neData = new SvNeoEpitope(fusion);
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

                // ignore IG and TR transcripts
                if(IMMUNE_TRANSCRIPT_PREFIXES.stream().anyMatch(x -> transData.BioType.startsWith(x)))
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
                    final String otherAminoAcidStr = otherNeData.aminoAcidString();

                    // remove exact matches or take the longer if one is a subset
                    if(aminoAcidStr.contains(otherAminoAcidStr))
                    {
                        neDataList.remove(j);
                    }
                    else if(otherAminoAcidStr.contains(aminoAcidStr))
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
            writeNeoepitopes(mWriter, mSampleId, false, neId, neData, upTransNames, downTransNames);
        }
    }
}
