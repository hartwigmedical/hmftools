package com.hartwig.hmftools.geneutils.mapping;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;

public class EnsemblTranscriptMapper
{
    private static final String ENSEMBL_DIR_SRC = "ensembl_dir_src";
    private static final String ENSEMBL_DIR_DEST = "ensembl_dir_dest";
    private static final String SPECIFIC_GENE_IDS = "specific_gene_ids";

    private final EnsemblDataCache mGeneCacheSrc;
    private final EnsemblDataCache mGeneCacheDest;
    private final List<String> mSpecificGeneIds;
    private final BufferedWriter mWriter;

    public EnsemblTranscriptMapper(final CommandLine cmd)
    {
        String outputDir = parseOutputDir(cmd);
        String ensemblDirSrc = cmd.getOptionValue(ENSEMBL_DIR_SRC);
        String ensemblDirDest = cmd.getOptionValue(ENSEMBL_DIR_DEST);

        if(outputDir == null || ensemblDirSrc == null || ensemblDirDest == null)
        {
            GU_LOGGER.error("missing config");
            System.exit(1);
        }

        mSpecificGeneIds = Lists.newArrayList();

        if(cmd.hasOption(SPECIFIC_GENE_IDS))
        {
            String specificGeneIds = cmd.getOptionValue(SPECIFIC_GENE_IDS);
            GU_LOGGER.info("loaded specific gene id(s): {}", specificGeneIds);
            Arrays.stream(specificGeneIds.split(";")).forEach(x -> mSpecificGeneIds.add(x));
        }

        // version makes no difference for transcript mapping since only appends 'chr' prefix
        mGeneCacheSrc = new EnsemblDataCache(ensemblDirSrc, RefGenomeVersion.V38);
        mGeneCacheSrc.setRequiredData(true, false, false, true);

        mGeneCacheDest = new EnsemblDataCache(ensemblDirDest, RefGenomeVersion.V38);
        mGeneCacheDest.setRequiredData(true, false, false, false);

        if(!mSpecificGeneIds.isEmpty())
        {
            mGeneCacheSrc.load(true);
            mGeneCacheSrc.loadTranscriptData(mSpecificGeneIds);

            mGeneCacheDest.load(true);
            mGeneCacheDest.loadTranscriptData(mSpecificGeneIds);
        }
        else
        {
            mGeneCacheSrc.load(false);
            mGeneCacheDest.load(false);
        }

        mGeneCacheSrc.createGeneIdDataMap();

        mWriter = initialiseWriter(outputDir);
    }

    public void run()
    {
        GU_LOGGER.info("running Ensembl transcript mapping");

        final Map<String, List<TranscriptData>> srcTrans = mGeneCacheSrc.getTranscriptDataMap();
        final Map<String, List<TranscriptData>> destTrans = mGeneCacheDest.getTranscriptDataMap();

        for(Map.Entry<String,List<TranscriptData>> geneEntry : mGeneCacheSrc.getTranscriptDataMap().entrySet())
        {
            String geneId = geneEntry.getKey();

            if(!mSpecificGeneIds.isEmpty() && !mSpecificGeneIds.contains(geneId))
                continue;

            GeneData geneData = mGeneCacheSrc.getGeneDataById(geneId);

            TranscriptData canonicalTransDataSrc = srcTrans.get(geneId).stream().filter(x -> x.IsCanonical).findFirst().orElse(null);

            if(canonicalTransDataSrc == null)
            {
                GU_LOGGER.error("gene({}: {}) canonical transcript not found", geneId, geneData.GeneName);
                continue;
            }

            if(canonicalTransDataSrc.CodingStart == null)
            {
                writeMappingData(geneData, canonicalTransDataSrc, Lists.newArrayList(), "NON_CODING");
                continue;
            }

            List<TranscriptData> transListDest = mGeneCacheDest.getTranscriptDataMap().get(geneId);

            if(transListDest == null || transListDest.isEmpty())
            {
                writeMappingData(geneData, canonicalTransDataSrc, Lists.newArrayList(), "UNMATCHED_GENE");
                continue;
            }
            else
            {
                mapTranscripts(geneData, canonicalTransDataSrc, transListDest);
            }
        }

        closeBufferedWriter(mWriter);
        GU_LOGGER.info("Ensembl transcript mapping complete");
    }

    private void mapTranscripts(final GeneData geneData, final TranscriptData transDataSrc, final List<TranscriptData> transDataDests)
    {
        // find any transcript which matches coding regions exactly
        int codingStart = transDataSrc.CodingStart;
        int codingEnd = transDataSrc.CodingEnd;

        final List<TranscriptData> matchedDestTrans = Lists.newArrayList();

        for(TranscriptData transDataDest : transDataDests)
        {
            if(transDataDest.CodingStart == null)
                continue;

            if(transDataDest.CodingStart != codingStart || transDataDest.CodingEnd != codingEnd)
                continue;

            // check exon boundaries within the coding region
            int exonDestIndex = 0;
            ExonData exonDest = null;

            while(exonDestIndex < transDataDest.exons().size())
            {
                exonDest = transDataDest.exons().get(exonDestIndex);

                if(exonDest.Start <= codingStart && exonDest.End >= codingStart)
                    break;

                ++exonDestIndex;
            }

            boolean matched = true;
            for(ExonData exonSrc : transDataSrc.exons())
            {
                if(exonSrc.End < codingStart)
                    continue;

                if(exonSrc.Start > codingEnd)
                    break;

                if(exonDestIndex >= transDataDest.exons().size())
                {
                    matched = false;
                    break;
                }

                exonDest = transDataDest.exons().get(exonDestIndex);

                int codingBaseStartSrc = max(codingStart, exonSrc.Start);
                int codingBaseEndSrc = min(codingEnd, exonSrc.End);

                int codingBaseStartDest = max(codingStart, exonDest.Start);
                int codingBaseEndDest = min(codingEnd, exonDest.End);

                if(codingBaseStartSrc != codingBaseStartDest || codingBaseEndSrc != codingBaseEndDest)
                {
                    matched = false;
                    break;
                }

                ++exonDestIndex;
            }

            if(matched)
                matchedDestTrans.add(transDataDest);
        }

        writeMappingData(geneData, transDataSrc, matchedDestTrans, !matchedDestTrans.isEmpty() ? "MATCHED" : "NO_MATCH");
    }

    private static BufferedWriter initialiseWriter(final String outputDir)
    {
        try
        {
            String outputFile = outputDir + "ensembl_transcript_mapping.csv";

            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("GeneId,GeneName,TransIdSrc,MatchInfo,MatchedCount,ContainsCanonical,ChangedCanonical,TransIdsDest");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            GU_LOGGER.error("error initialising Ensembl transcript mapping file: {}", e.toString());
            return null;
        }
    }

    private void writeMappingData(
            final GeneData geneData, final TranscriptData transDataSrc,
            final List<TranscriptData> transDataDests, final String matchInfo)
    {
        try
        {
            mWriter.write(String.format("%s,%s,%s,%s",
                    geneData.GeneId, geneData.GeneName, transDataSrc.TransName, matchInfo));

            if(transDataDests.isEmpty())
            {
                mWriter.write(",0,false,true,");
            }
            else
            {
                TranscriptData newCanonical = transDataDests.stream().filter(x -> x.IsCanonical).findFirst().orElse(null);
                boolean changedCanonical = newCanonical == null || !newCanonical.TransName.equals(transDataSrc.TransName);
                StringJoiner sj = new StringJoiner(";");
                transDataDests.forEach(x -> sj.add(x.TransName));

                mWriter.write(String.format(",%s,%s,%s,%s",
                        transDataDests.size(), newCanonical != null, changedCanonical, sj.toString()));
            }

            mWriter.newLine();
        }
        catch(IOException e)
        {
            GU_LOGGER.error("error writing Ensembl transcript mapping file: {}", e.toString());
        }
    }

    public static void main(String[] args) throws ParseException
    {
        final Options options = createOptions();
        final CommandLine cmd = new DefaultParser().parse(options, args);

        setLogLevel(cmd);

        EnsemblTranscriptMapper transcriptMapper = new EnsemblTranscriptMapper(cmd);
        transcriptMapper.run();
    }

    private static Options createOptions()
    {
        final Options options = new Options();
        options.addOption(RefGenomeVersion.REF_GENOME_VERSION, true, "Ref genome version (V37 or V38))");
        options.addOption(ENSEMBL_DIR_SRC, true, "Ensembl data cache dir for ref-genome v37");
        options.addOption(ENSEMBL_DIR_DEST, true, "Ensembl data cache dir for ref-genome v38");
        options.addOption(SPECIFIC_GENE_IDS, true, "Optional list of geneIds to map");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(LOG_DEBUG, false, "Log verbose");
        return options;
    }

}
