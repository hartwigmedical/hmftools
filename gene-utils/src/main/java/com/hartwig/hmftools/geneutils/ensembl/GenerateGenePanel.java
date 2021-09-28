package com.hartwig.hmftools.geneutils.ensembl;

import static com.hartwig.hmftools.common.genome.genepanel.HmfTranscriptRegionFile.DEFAULT_DELIM;
import static com.hartwig.hmftools.common.genome.genepanel.HmfTranscriptRegionFile.EXON_DATA_DELIM;
import static com.hartwig.hmftools.common.genome.genepanel.HmfTranscriptRegionFile.ITEM_DELIM;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.genepanel.HmfTranscriptRegionFile;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class GenerateGenePanel
{
    private static final String CDKN2A_ALT = "ENST00000579755";

    // could consider loading these from the alt-transcripts soon to be in the driver panel file
    private static final List<String> NON_CANONICAL_TRANSCRIPTS = Lists.newArrayList(CDKN2A_ALT);

    public static void main(String[] args) throws ParseException, IOException, SQLException
    {
        final Options options = createOptions();
        final CommandLine cmd = new DefaultParser().parse(options, args);

        setLogLevel(cmd);

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.from(cmd.getOptionValue(RefGenomeVersion.REF_GENOME_VERSION));
        String outputDir = parseOutputDir(cmd);

        GU_LOGGER.info("writing Ensembl gene panel, ref-genome-version({})", refGenomeVersion);

        String geneFile = String.format("%sall_genes.%s.tsv", outputDir, refGenomeVersion.identifier());

        EnsemblDataCache ensemblDataCache = new EnsemblDataCache(cmd, refGenomeVersion);
        ensemblDataCache.setRequiredData(true, false, false, false);
        ensemblDataCache.load(false);

        try
        {
            BufferedWriter writer = createBufferedWriter(geneFile, false);

            writer.write(HmfTranscriptRegionFile.header());
            writer.newLine();

            for(Map.Entry<String, List<GeneData>> entry : ensemblDataCache.getChrGeneDataMap().entrySet())
            {
                for(final GeneData geneData : entry.getValue())
                {
                    List<TranscriptData> transDataList = ensemblDataCache.getTranscripts(geneData.GeneId);

                    if(transDataList == null)
                        continue;

                    for(TranscriptData transData : transDataList)
                    {
                        if(!transData.IsCanonical && !NON_CANONICAL_TRANSCRIPTS.contains(transData.TransName))
                            continue;

                        StringJoiner sj = new StringJoiner(DEFAULT_DELIM);
                        sj.add(geneData.GeneId);
                        sj.add(geneData.GeneName);
                        sj.add(refGenomeVersion.versionedChromosome(geneData.Chromosome));
                        sj.add(String.valueOf(geneData.GeneStart));
                        sj.add(String.valueOf(geneData.GeneEnd));
                        sj.add(String.valueOf(geneData.Strand));
                        sj.add(geneData.getSynonyms());
                        sj.add(geneData.KaryotypeBand);

                        sj.add(transData.TransName);
                        sj.add(String.valueOf(transData.IsCanonical));
                        sj.add(String.valueOf(transData.TransStart));
                        sj.add(String.valueOf(transData.TransEnd));
                        sj.add(transData.CodingStart != null ? String.valueOf(transData.CodingStart) : "");
                        sj.add(transData.CodingEnd != null ? String.valueOf(transData.CodingEnd) : "");

                        StringJoiner exonSj = new StringJoiner(ITEM_DELIM);
                        transData.exons().forEach(x -> exonSj.add(String.format("%d%s%d", x.Start, EXON_DATA_DELIM, x.End)));
                        sj.add(exonSj.toString());

                        writer.write(sj.toString());
                        writer.newLine();
                    }
                }
            }

            writer.close();

            GU_LOGGER.info("written gene panel output to {}", geneFile);
        }
        catch (final IOException e)
        {
            GU_LOGGER.error("error writing gene panel data file: {}", e.toString());
            return;
        }
    }

    private static Options createOptions()
    {
        final Options options = new Options();
        options.addOption(RefGenomeVersion.REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);
        EnsemblDataCache.addEnsemblDir(options);
        addLoggingOptions(options);
        addOutputDir(options);
        return options;
    }

}
