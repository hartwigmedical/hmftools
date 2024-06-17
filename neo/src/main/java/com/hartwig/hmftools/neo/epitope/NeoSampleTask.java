package com.hartwig.hmftools.neo.epitope;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.AminoAcidRna.STOP_SYMBOL;
import static com.hartwig.hmftools.common.gene.TranscriptProteinData.BIOTYPE_NONSENSE_MED_DECAY;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFusion.generateFilename;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact.VAR_TRANS_IMPACT_ANNOTATION;
import static com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact.fromVariantContext;
import static com.hartwig.hmftools.neo.NeoCommon.DOWNSTREAM_PRE_GENE_DISTANCE;
import static com.hartwig.hmftools.neo.NeoCommon.IMMUNE_TRANSCRIPT_PREFIXES;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.epitope.NeoEpitopeFinder.initialiseNeoepitopeWriter;
import static com.hartwig.hmftools.neo.epitope.NeoEpitopeFinder.writeNeoepitopes;
import static com.hartwig.hmftools.neo.epitope.PointMutationData.checkVariantEffects;
import static com.hartwig.hmftools.neo.epitope.PointMutationData.isRelevantMutation;
import static com.hartwig.hmftools.neo.epitope.SvNeoEpitope.svIsNonDisruptiveInCodingTranscript;

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
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;

import htsjdk.variant.variantcontext.VariantContext;

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

        String somaticVcf = convertWildcardSamplePath(mConfig.SomaticVcf, mSampleId);

        if(!Files.exists(Paths.get(somaticVcf)))
        {
            NE_LOGGER.warn("Purple somatic VCF file({}) not found", somaticVcf);
            return pointMutations;
        }

        VcfFileReader reader = new VcfFileReader(somaticVcf);

        if(!reader.fileValid())
        {
            System.exit(1);
        }

        try
        {
            for(VariantContext variantContext : reader.iterator())
            {
                if(variantContext.isFiltered())
                    continue;

                VariantContextDecorator variant = new VariantContextDecorator(variantContext);

                // must have a gene impact
                boolean hasValidGeneImpact = false;
                String geneName = "";
                CodingEffect worstCodingEffect = NONE;

                VariantImpact variantImpact = variant.variantImpact();

                if(!variantImpact.GeneName.isEmpty())
                {
                    if(mGeneTransCache.getGeneDataByName(variantImpact.GeneName) != null && isRelevantMutation(variantImpact))
                    {
                        worstCodingEffect = variantImpact.WorstCodingEffect;
                        geneName = variantImpact.GeneName;
                        hasValidGeneImpact = true;
                    }
                }

                if(!hasValidGeneImpact && variantContext.hasAttribute(VAR_TRANS_IMPACT_ANNOTATION))
                {
                    List<VariantTranscriptImpact> transImpacts = fromVariantContext(variantContext);

                    for(VariantTranscriptImpact transcriptImpact : transImpacts)
                    {
                        if(mGeneTransCache.getGeneDataByName(transcriptImpact.GeneName) == null)
                            continue;

                        CodingEffect codingEffect = checkVariantEffects(transcriptImpact);

                        if(codingEffect != NONE)
                        {
                            worstCodingEffect = codingEffect;
                            geneName = transcriptImpact.GeneName;
                            hasValidGeneImpact = true;
                            break;
                        }
                    }
                }

                if(!hasValidGeneImpact)
                    continue;

                pointMutations.add(new PointMutationData(
                        variant.chromosome(), variant.position(), variant.ref(), variant.alt(), geneName, worstCodingEffect,
                        variant.variantCopyNumber(), variant.adjustedCopyNumber(),
                        variantContext.getAttributeAsDouble(SUBCLONAL_LIKELIHOOD_FLAG, 0)));
            }

            NE_LOGGER.debug("loaded {} somatic variants from file({})", pointMutations.size(), somaticVcf);
        }
        catch(Exception e)
        {
            NE_LOGGER.error(" failed to read somatic VCF file({}): {}", somaticVcf, e.toString());
            System.exit(1);
        }

        return pointMutations;
    }

    private List<NeoEpitopeFusion> getSvFusions()
    {
        List<NeoEpitopeFusion> fusions = Lists.newArrayList();

        if(mConfig.LinxDir == null)
            return fusions;

        String linxDir = convertWildcardSamplePath(mConfig.LinxDir, mSampleId);
        final String filename = generateFilename(linxDir, mSampleId);

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

            // for same gene fusions, check that the SV isn't non-disruptive within the coding region of any other transcript(s)

            boolean isNonDisruptiveCoding = sameGene
                    && fusion.Orientations[FS_UP] == POS_ORIENT && fusion.Orientations[FS_DOWN] == NEG_ORIENT
                    && upTransDataList.stream().anyMatch(x -> svIsNonDisruptiveInCodingTranscript(fusion.Positions, x));

            if(isNonDisruptiveCoding)
                continue;

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
