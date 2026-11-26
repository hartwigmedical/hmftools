package com.hartwig.hmftools.qsee.coverage;

import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.qsee.common.QseeFileCommon;

import org.jetbrains.annotations.Nullable;

public class ExonExclusions
{
    public static final String EXON_EXCLUSIONS_FILE_CFG = "exon_exclusions_file";
    public static final String EXON_EXCLUSIONS_FILE_CFG_DESC = "Path to exon exclusions file with columns GeneName and ExonRank (one row per exon)";

    private static final String COL_GENE_NAME = "GeneName";
    private static final String COL_EXON_RANK = "ExonRank";

    private final Map<String,Set<Integer>> mExonExclusions;

    public ExonExclusions(Map<String,Set<Integer>> exonExclusions)
    {
        mExonExclusions = exonExclusions;
    }

    public ExonExclusions()
    {
        mExonExclusions = new HashMap<>();
    }

    public void add(String geneName, int exonRank)
    {
        mExonExclusions.computeIfAbsent(geneName, k -> new HashSet<>()).add(exonRank);
    }

    public boolean isExcluded(String geneName, int exonRank)
    {
        return mExonExclusions.getOrDefault(geneName, Set.of()).contains(exonRank);
    }

    public static ExonExclusions empty() { return new ExonExclusions(); }

    public static ExonExclusions read(String filename)
    {
        QC_LOGGER.info("Reading exon exclusions from: {}", filename);

        ExonExclusions exonExclusions = new ExonExclusions();

        try(DelimFileReader reader = new DelimFileReader(filename))
        {
            for(DelimFileReader.Row row : reader)
            {
                String geneName = row.get(COL_GENE_NAME);
                int exonRank = row.getInt(COL_EXON_RANK);
                exonExclusions.add(geneName, exonRank);
            }
        }

        return exonExclusions;
    }

    public static String generateFilename(String outputDir, @Nullable String outputId)
    {
        return QseeFileCommon.generateFilename(outputDir, null, "exon_exclusions", outputId, "tsv");
    }

    public void write(String filename)
    {
        QC_LOGGER.info("Writing exon exclusions to: {}", filename);

        List<Map.Entry<String,Integer>> excludedExons = new ArrayList<>();
        for(String geneName : new TreeSet<>(mExonExclusions.keySet()))
        {
            for(int exonRank : new TreeSet<>(mExonExclusions.get(geneName)))
                excludedExons.add(Map.entry(geneName, exonRank));
        }

        DelimFileWriter.write(filename, List.of(COL_GENE_NAME, COL_EXON_RANK), excludedExons, (excludedExon, row) -> {
            row.set(COL_GENE_NAME, excludedExon.getKey());
            row.set(COL_EXON_RANK, excludedExon.getValue());
        });
    }
}
