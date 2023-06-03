package com.hartwig.hmftools.geneutils.mapping;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;

public class EnsemblCacheCompare
{
    private static final String ENSEMBL_DIR_REF = "ensembl_dir_ref";
    private static final String ENSEMBL_DIR_NEW = "ensembl_dir_new";
    private static final String SPECIFIC_GENES = "specific_genes";

    private final EnsemblDataCache mGeneCacheRef;
    private final EnsemblDataCache mGeneCacheNew;
    private final List<String> mSpecificGeneNames;
    private final BufferedWriter mWriter;

    public EnsemblCacheCompare(final CommandLine cmd)
    {
        String outputDir = parseOutputDir(cmd);
        String outputId = cmd.getOptionValue(OUTPUT_ID);
        String ensemblDirRef = cmd.getOptionValue(ENSEMBL_DIR_REF);
        String ensemblDirNew = cmd.getOptionValue(ENSEMBL_DIR_NEW);

        if(outputDir == null || ensemblDirRef == null || ensemblDirNew == null)
        {
            GU_LOGGER.error("missing config");
            System.exit(1);
        }

        mSpecificGeneNames = Lists.newArrayList();

        if(cmd.hasOption(SPECIFIC_GENES))
        {
            String specificGeneIds = cmd.getOptionValue(SPECIFIC_GENES);
            GU_LOGGER.info("loaded specific gene id(s): {}", specificGeneIds);
            Arrays.stream(specificGeneIds.split(";")).forEach(x -> mSpecificGeneNames.add(x));
        }

        // version makes no difference for transcript mapping since only appends 'chr' prefix
        mGeneCacheRef = new EnsemblDataCache(ensemblDirRef, RefGenomeVersion.V37);
        mGeneCacheRef.setRequiredData(true, false, false, false);

        mGeneCacheNew = new EnsemblDataCache(ensemblDirNew, RefGenomeVersion.V37);
        mGeneCacheNew.setRequiredData(true, false, false, false);

        if(!mSpecificGeneNames.isEmpty())
        {
            mGeneCacheRef.load(true);
            mGeneCacheRef.loadTranscriptData(mSpecificGeneNames);

            mGeneCacheNew.load(true);
            mGeneCacheNew.loadTranscriptData(mSpecificGeneNames);
        }
        else
        {
            mGeneCacheRef.load(false);
            mGeneCacheNew.load(false);
        }

        mGeneCacheRef.createGeneIdDataMap();
        mGeneCacheNew.createGeneIdDataMap();

        mWriter = initialiseWriter(outputDir, outputId);
    }

    private enum DiffScope
    {
        GENE,
        TRANSCRIPT,
        EXON;
    }

    private enum MismatchType
    {
        REF_ONLY,
        NEW_ONLY,
        VALUE;
    }

    private enum TranscriptDiff
    {
        CODING_NON_CODING,
        CODING_REGION,
        BIOTYPE,
        EXONS;
    }

    public void run()
    {
        GU_LOGGER.info("running Ensembl cache comparisons");

        // types of differences:
        // gene differences - new genes added, others no longer present
        // with each gene - transcript differences

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = chromosome.toString();
            List<GeneData> geneDataListRef = mGeneCacheRef.getChrGeneDataMap().get(chrStr);
            List<GeneData> geneDataListNew = mGeneCacheNew.getChrGeneDataMap().get(chrStr);

            if(geneDataListNew == null || geneDataListRef == null)
            {
                GU_LOGGER.error("chromosome({}) missing gene data: refValid({}) newValid({})",
                        chrStr, geneDataListRef != null, geneDataListNew != null);
                continue;
            }

            GU_LOGGER.debug("chromosome({}) comparing gene data, counts ref({}) new({})",
                    chrStr, geneDataListRef.size(), geneDataListNew.size());

            Set<String> newProcessed = Sets.newHashSet();
            for(GeneData geneDataRef : geneDataListRef)
            {
                GeneData geneDataNew = mGeneCacheNew.getGeneDataByName(geneDataRef.GeneName);

                compareGenes(geneDataRef, geneDataNew);
                newProcessed.add(geneDataRef.GeneName);
            }

            for(GeneData geneDataNew : geneDataListNew)
            {
                if(newProcessed.contains(geneDataNew.GeneName))
                    continue;

                compareGenes(null, geneDataNew);
            }
        }

        closeBufferedWriter(mWriter);
        GU_LOGGER.info("Ensembl cache comparison complete");
    }

    private void compareGenes(final GeneData geneDataRef, final GeneData geneDataNew)
    {
        if(!mSpecificGeneNames.isEmpty())
        {
            String geneName = geneDataRef != null ? geneDataRef.GeneName : geneDataNew.GeneName;

            if(!mSpecificGeneNames.contains(geneName))
                return;
        }

        if(geneDataNew == null || geneDataRef == null)
        {
            writeGeneMismatch(geneDataRef, geneDataNew, "");
            return;
        }

        // compare transcripts
        List<TranscriptData> transcriptsRef = mGeneCacheRef.getTranscriptDataMap().get(geneDataRef.GeneId);
        List<TranscriptData> transcriptsNew = mGeneCacheNew.getTranscriptDataMap().get(geneDataNew.GeneId);

        Set<String> processedNew = Sets.newHashSet();

        if(transcriptsRef != null)
        {
            for(TranscriptData transDataRef : transcriptsRef)
            {
                TranscriptData transDataNew = transcriptsNew != null ?
                        transcriptsNew.stream().filter(x -> x.TransName.equals(transDataRef.TransName)).findFirst().orElse(null) : null;

                if(transDataNew == null)
                {
                    writeTranscriptMismatch(geneDataRef, transDataRef, transDataNew, "");
                    continue;
                }
                else
                {
                    processedNew.add(transDataNew.TransName);

                    // compare transcript details
                    compareTranscripts(geneDataRef, transDataRef, transDataNew);
                }
            }
        }

        if(transcriptsNew != null)
        {
            for(TranscriptData transDataNew : transcriptsNew)
            {
                if(processedNew.contains(transDataNew.TransName))
                    continue;

                writeTranscriptMismatch(geneDataRef, null, transDataNew, "");
            }
        }
    }

    private void compareTranscripts(final GeneData geneData, final TranscriptData transDataRef, final TranscriptData transDataNew)
    {
        List<TranscriptDiff> diffs = Lists.newArrayList();

        if(transDataRef.nonCoding() != transDataNew.nonCoding())
        {
            diffs.add(TranscriptDiff.CODING_NON_CODING);
        }

        if(transDataRef.exons().size() != transDataNew.exons().size())
        {
            diffs.add(TranscriptDiff.EXONS);
        }

        if(!transDataRef.BioType.equals(transDataNew.BioType))
        {
            diffs.add(TranscriptDiff.BIOTYPE);
        }

        if(!transDataNew.nonCoding() && !transDataRef.nonCoding())
        {
            if(!transDataRef.CodingStart.equals(transDataNew.CodingStart) || !transDataRef.CodingEnd.equals(transDataNew.CodingEnd))
            {
                diffs.add(TranscriptDiff.CODING_REGION);
            }
        }

        if(!diffs.isEmpty())
        {
            StringJoiner sj = new StringJoiner(";");
            diffs.forEach(x -> sj.add(x.toString()));
            writeTranscriptMismatch(geneData, transDataRef, transDataNew, sj.toString());
        }
    }

    private static BufferedWriter initialiseWriter(final String outputDir, final String outputId)
    {
        try
        {
            String outputFile = outputDir + "ensembl_cache_diffs." + outputId + TSV_EXTENSION;

            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("Scope\tGeneName\tGeneId\tTransName\tMismatchType\tInfo\tDiffs");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            GU_LOGGER.error("error initialising Ensembl cache comparison file: {}", e.toString());
            return null;
        }
    }

    private void writeGeneMismatch(
            final GeneData geneDataRef, final GeneData geneDataNew, final String diffs)
    {
        try
        {
            GeneData validGeneData = geneDataRef != null ? geneDataRef : geneDataNew;
            MismatchType mismatchType = geneDataNew == null ? MismatchType.REF_ONLY : MismatchType.NEW_ONLY;

            String geneInfo = format("loc(%s:%d-%d) strand(%d)",
                    validGeneData.Chromosome, validGeneData.GeneStart, validGeneData.GeneEnd, validGeneData.Strand);

            mWriter.write(format("%s\t%s\t%s\t%s\t%s\t%s\t%s",
                    DiffScope.GENE, validGeneData.GeneName, validGeneData.GeneId, "", mismatchType, geneInfo, diffs));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            GU_LOGGER.error("error writing Ensembl cache comparison file: {}", e.toString());
        }
    }

    private void writeTranscriptMismatch(
            final GeneData geneData, final TranscriptData transDataRef, final TranscriptData transDataNew, final String diffs)
    {
        try
        {
            TranscriptData validTransData = transDataRef != null ? transDataRef : transDataNew;
            MismatchType mismatchType;

            if(transDataRef != null && transDataNew != null)
                mismatchType = MismatchType.VALUE;
            else if(transDataNew == null)
                mismatchType = MismatchType.REF_ONLY;
            else
                mismatchType = MismatchType.NEW_ONLY;

            String transInfo = format("loc(%s:%d-%d) strand(%d) protein(%s/%s) exons(%d/%d) codingRegion(%d:%d/%d:%d)",
                    geneData.Chromosome, validTransData.TransStart, validTransData.TransEnd, validTransData.Strand,
                    transDataRef != null ? transDataRef.BioType : "", transDataNew != null ? transDataNew.BioType : "",
                    transDataRef != null ? transDataRef.exons().size() : 0, transDataNew != null ? transDataNew.exons().size() : 0,
                    transDataRef != null &&  transDataRef.CodingStart != null ? transDataRef.CodingStart : 0,
                    transDataRef != null &&  transDataRef.CodingEnd != null ? transDataRef.CodingEnd : 0,
                    transDataNew != null &&  transDataNew.CodingStart != null ? transDataNew.CodingStart : 0,
                    transDataNew != null &&  transDataNew.CodingEnd != null ? transDataNew.CodingEnd : 0);

            mWriter.write(format("%s\t%s\t%s\t%s\t%s\t%s\t%s",
                    DiffScope.TRANSCRIPT, geneData.GeneName, validTransData.GeneId, validTransData.TransName, mismatchType, transInfo, diffs));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            GU_LOGGER.error("error writing Ensembl cache comparison file: {}", e.toString());
        }
    }

    private void writeMappingData(
            final GeneData geneData, final TranscriptData transDataSrc,
            final List<TranscriptData> transDataDests, final String matchInfo)
    {
        try
        {
            mWriter.write(format("%s,%s,%s,%s",
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

                mWriter.write(format(",%s,%s,%s,%s",
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

        EnsemblCacheCompare cacheCompare = new EnsemblCacheCompare(cmd);
        cacheCompare.run();
    }

    private static Options createOptions()
    {
        final Options options = new Options();
        // options.addOption(RefGenomeVersion.REF_GENOME_VERSION, true, "Ref genome version (V37 or V38))");
        options.addOption(ENSEMBL_DIR_REF, true, "Ensembl data cache dir for ref-genome v37");
        options.addOption(ENSEMBL_DIR_NEW, true, "Ensembl data cache dir for ref-genome v38");
        options.addOption(SPECIFIC_GENES, true, "Limit comparison to specific genes, separated by ','");
        addOutputOptions(options);
        addLoggingOptions(options);
        return options;
    }

}
