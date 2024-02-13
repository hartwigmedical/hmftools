package com.hartwig.hmftools.geneutils.ensembl;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.ENSEMBL_TRANS_AMINO_ACIDS_FILE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.calcCodingStartPositionAdjustment;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import org.jetbrains.annotations.NotNull;

public class ProteomeWriter
{
    private final EnsemblDataCache mEnsemblDataCache;
    private final RefGenomeInterface mRefGenome;
    private final RefGenomeVersion mRefGenomeVersion;
    private final boolean mCanonicalOnly;

    private static final String CANONICAL_ONLY = "canonical_only";
    private static final String GENE_IDS = "gene_ids";

    private BufferedWriter mWriter;

    public ProteomeWriter(final ConfigBuilder configBuilder)
    {
        mCanonicalOnly = configBuilder.hasFlag(CANONICAL_ONLY);

        mRefGenomeVersion = RefGenomeVersion.from(configBuilder);
        mEnsemblDataCache = new EnsemblDataCache(configBuilder);
        mEnsemblDataCache.setRequiredData(true, false, false, mCanonicalOnly);

        if(configBuilder.hasValue(GENE_IDS))
        {
            List<String> specificGeneIds = Lists.newArrayList();
            Arrays.stream(configBuilder.getValue(GENE_IDS).split(",", -1)).forEach(x -> specificGeneIds.add(x));
            mEnsemblDataCache.setRestrictedGeneIdList(specificGeneIds);
        }

        mRefGenome = loadRefGenome(configBuilder.getValue(REF_GENOME));

        mWriter = initialiseWriter(parseOutputDir(configBuilder));
    }

    public void run()
    {
        mEnsemblDataCache.load(false);

        GU_LOGGER.info("writing proteome to file");

        int geneCount = 0;

        for(Map.Entry<String,List<GeneData>> entry : mEnsemblDataCache.getChrGeneDataMap().entrySet())
        {
            for(GeneData geneData : entry.getValue())
            {
                List<TranscriptData> transDataList = mEnsemblDataCache.getTranscripts(geneData.GeneId);

                if(transDataList == null)
                    continue;

                for(TranscriptData transData : transDataList)
                {
                    if(mCanonicalOnly && !transData.IsCanonical)
                        continue;

                    processTranscript(geneData, transData);

                    if(mCanonicalOnly)
                        break;
                }

                ++geneCount;
            }
        }

        closeBufferedWriter(mWriter);

        GU_LOGGER.info("wrote {} gene amino-acids sequences", geneCount);

        GU_LOGGER.info("proteome write complete");
    }

    private void processTranscript(final GeneData geneData, final TranscriptData transData)
    {
        if(transData.CodingStart == null)
            return;

        boolean inCoding = false;
        String aminoAcids = "";

        String chromosome = mRefGenomeVersion.versionedChromosome(geneData.Chromosome);

        if(transData.Strand == POS_STRAND)
        {
            StringBuilder codingBases = new StringBuilder();

            for(int i = 0; i < transData.exons().size(); ++i)
            {
                ExonData exon = transData.exons().get(i);

                if(exon.End < transData.CodingStart)
                    continue;

                if(exon.Start > transData.CodingEnd)
                    break;

                int exonCodingStart = max(transData.CodingStart, exon.Start);
                int exonCodingEnd = min(transData.CodingEnd, exon.End);

                if(!inCoding)
                {
                    inCoding = true;
                    exonCodingStart += calcCodingStartPositionAdjustment(transData, exon);
                }

                codingBases.append(mRefGenome.getBaseString(chromosome, exonCodingStart, exonCodingEnd));
            }

            aminoAcids = Codons.aminoAcidFromBases(codingBases.toString());
        }
        else
        {
            String codingBases = "";
            for(int i = transData.exons().size() - 1; i >= 0; --i)
            {
                ExonData exon = transData.exons().get(i);

                if(exon.Start > transData.CodingEnd)
                    continue;

                if(exon.End < transData.CodingStart)
                    break;

                int exonCodingStart = max(transData.CodingStart, exon.Start);
                int exonCodingEnd = min(transData.CodingEnd, exon.End);

                if(!inCoding)
                {
                    inCoding = true;
                    exonCodingEnd += calcCodingStartPositionAdjustment(transData, exon);
                }

                codingBases = mRefGenome.getBaseString(chromosome, exonCodingStart, exonCodingEnd) + codingBases;
            }

            aminoAcids = Codons.aminoAcidFromBases(Nucleotides.reverseComplementBases(codingBases));
        }

        writeData(geneData.GeneId, geneData.GeneName, transData.TransName, transData.IsCanonical, aminoAcids);
    }

    private static BufferedWriter initialiseWriter(final String outputDir)
    {
        try
        {
            String filename = outputDir + ENSEMBL_TRANS_AMINO_ACIDS_FILE;
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write(TranscriptAminoAcids.csvHeader());
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to initialise CSV file output: {}", e.toString());
            return null;
        }
    }

    private void writeData(final String geneId, final String geneName, final String transName, boolean isCanonical, final String aminoAcids)
    {
        try
        {
            mWriter.write(String.format("%s,%s,%s,%s,%s", geneId, geneName, transName, isCanonical, aminoAcids));
            mWriter.newLine();
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write peptide data: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        addEnsemblDir(configBuilder, true);
        addRefGenomeConfig(configBuilder, true);
        addOutputDir(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
        configBuilder.addFlag(CANONICAL_ONLY, "Only write canonical proteome");
        configBuilder.addConfigItem(GENE_IDS, "Specific gene IDs only");

        configBuilder.checkAndParseCommandLine(args);

        ProteomeWriter proteomeWriter = new ProteomeWriter(configBuilder);
        proteomeWriter.run();
    }
}
