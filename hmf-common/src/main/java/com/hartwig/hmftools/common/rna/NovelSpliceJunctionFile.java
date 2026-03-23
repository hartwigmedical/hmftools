package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_DEPTH_END;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_DEPTH_START;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_FRAG_COUNT;
import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.RNA_LOGGER;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_ID;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.Nullable;

public final class NovelSpliceJunctionFile
{
    public static final String ALT_SJ_FILE_ID = "alt_splice_junc.tsv";
    public static final String ALT_SJ_UNFILTERED_FILE_ID = "alt_splice_junc_unfiltered.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ISF_FILE_ID + ALT_SJ_FILE_ID;
    }

    public static String generateUnfilteredFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ISF_FILE_ID + ALT_SJ_UNFILTERED_FILE_ID;
    }

    private enum Columns
    {
        GeneName,
        Chromosome,
        SjStart,
        SjEnd,
        Type,
        TranscriptStart,
        TranscriptEnd,
        ExonStart,
        ExonEnd,
        FragCount,
        DepthStart,
        DepthEnd,
        RegionStart,
        RegionEnd,
        BaseStart,
        BaseEnd,
        CohortFrequency;
    }

    public static final String FLD_ALT_SJ_POS_START = Columns.SjStart.toString();
    public static final String FLD_ALT_SJ_POS_END = Columns.SjEnd.toString();
    public static final String FLD_ALT_SJ_TYPE = Columns.Type.toString();
    public static final String FLD_BASES_START = Columns.BaseStart.toString();
    public static final String FLD_BASES_END = Columns.BaseEnd.toString();
    public static final String FLD_COHORT_FREQUENCY = Columns.CohortFrequency.toString();

    public static String formKey(final String chromosome, int posStart, int posEnd)
    {
        return String.format("%s;%s;%s", chromosome, posStart, posEnd);
    }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        Arrays.stream(Columns.values()).forEach(x -> sj.add(x.toString()));
        return sj.toString();
    }

    public static String write(final NovelSpliceJunction novelSpliceJunction)
    {
        return new StringJoiner(TSV_DELIM)
                .add(novelSpliceJunction.geneName())
                .add(novelSpliceJunction.chromosome())
                .add(String.valueOf(novelSpliceJunction.junctionStart()))
                .add(String.valueOf(novelSpliceJunction.junctionEnd()))
                .add(novelSpliceJunction.type().toString())
                .add(novelSpliceJunction.transcriptStart())
                .add(novelSpliceJunction.transcriptEnd())
                .add(String.valueOf(novelSpliceJunction.exonStart()))
                .add(String.valueOf(novelSpliceJunction.exonEnd()))
                .add(String.valueOf(novelSpliceJunction.fragmentCount()))
                .add(String.valueOf(novelSpliceJunction.depthStart()))
                .add(String.valueOf(novelSpliceJunction.depthEnd()))
                .add(novelSpliceJunction.regionStart().toString())
                .add(novelSpliceJunction.regionEnd().toString())
                .add(novelSpliceJunction.basesStart())
                .add(novelSpliceJunction.basesEnd())
                .add(String.valueOf(novelSpliceJunction.cohortFrequency()))
                .toString();
    }

    public static List<NovelSpliceJunction> read(final String filename)
    {
        try
        {
            List<NovelSpliceJunction> novelJunctions = Lists.newArrayList();

            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();
            String fileDelim = inferFileDelimiter(filename);
            final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(line, fileDelim);

            while((line = fileReader.readLine()) != null)
            {
                String[] values = line.split(fileDelim, -1);

                int geneIndex = fieldsIndexMap.get(Columns.GeneName.toString());
                int chrIndex = fieldsIndexMap.get(Columns.Chromosome.toString());
                int sjStartIndex = fieldsIndexMap.get(Columns.SjStart.toString());
                int sjEndIndex = fieldsIndexMap.get(Columns.SjEnd.toString());
                int typeIndex = fieldsIndexMap.get(Columns.Type.toString());
                int fragCountIndex = fieldsIndexMap.get(Columns.FragCount.toString());
                int depthStartIndex = fieldsIndexMap.get(Columns.DepthStart.toString());
                int depthEndIndex = fieldsIndexMap.get(Columns.DepthEnd.toString());
                int regionStartIndex = fieldsIndexMap.get(Columns.RegionStart.toString());
                int regionEndIndex = fieldsIndexMap.get(Columns.RegionEnd.toString());
                Integer basesStartIndex = fieldsIndexMap.get(Columns.BaseStart.toString());
                Integer basesEndIndex = fieldsIndexMap.get(Columns.BaseEnd.toString());
                Integer cohortFreqIndex = fieldsIndexMap.get(Columns.CohortFrequency.toString());

                Integer transStartIndex = fieldsIndexMap.get(Columns.TranscriptStart.toString());
                Integer transEndIndex = fieldsIndexMap.get(Columns.TranscriptEnd.toString());
                Integer exonStartIndex = fieldsIndexMap.get(Columns.ExonStart.toString());
                Integer exonEndIndex = fieldsIndexMap.get(Columns.ExonEnd.toString());

                String chromosome = values[chrIndex];
                int sjStart = Integer.parseInt(values[sjStartIndex]);
                int sjEnd = Integer.parseInt(values[sjEndIndex]);

                novelJunctions.add(ImmutableNovelSpliceJunction.builder()
                        .geneName(values[geneIndex])
                        .chromosome(chromosome)
                        .junctionStart(sjStart)
                        .junctionEnd(sjEnd)
                        .type(AltSpliceJunctionType.valueOf(values[typeIndex]))
                        .transcriptStart(transStartIndex != null ? values[transStartIndex] : "")
                        .transcriptEnd(transEndIndex != null ? values[transEndIndex] : "")
                        .exonStart(exonStartIndex != null ? Integer.parseInt(values[exonStartIndex]) : -1)
                        .exonEnd(exonEndIndex != null ? Integer.parseInt(values[exonEndIndex]) : -1)
                        .fragmentCount(Integer.parseInt(values[fragCountIndex]))
                        .depthStart(Integer.parseInt(values[depthStartIndex]))
                        .depthEnd(Integer.parseInt(values[depthEndIndex]))
                        .regionStart(AltSpliceJunctionContext.valueOf(values[regionStartIndex]))
                        .regionEnd(AltSpliceJunctionContext.valueOf(values[regionEndIndex]))
                        .basesStart(basesStartIndex != null ? values[basesStartIndex] : "")
                        .basesEnd(basesEndIndex != null ? values[basesEndIndex] : "")
                        .cohortFrequency(cohortFreqIndex != null ? Integer.parseInt(values[cohortFreqIndex]) : 0)
                        .build());
            }

            return novelJunctions;
        }
        catch(IOException e)
        {
            RNA_LOGGER.error("failed to load Isofox fusion file({}): {}", filename, e.toString());
            return null;
        }
    }
}
